
import logging
import os
import re
import shutil  #change to use xml from blast web
import subprocess
import time
from io import StringIO
from multiprocessing import cpu_count
from tempfile import NamedTemporaryFile
from textwrap import dedent

import requests
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

from .. import parameters_default
from .obtaining_sequences import get_seq


logger = logging.getLogger()

# Remote BLAST (NCBI) settings
NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
NCBI_TOOL = 'ProLink'
NCBI_EMAIL = 'example@example.com'  # contact email sent to NCBI (see the NCBI usage guidelines)
BLAST_POLL_INTERVAL = 60  # seconds between status checks (NCBI asks not to poll a single RID more than once a minute)
BLAST_TIMEOUT = 3600      # seconds to wait for a single search before giving up (def: 1 hour)
BLAST_MAX_ATTEMPTS = 3    # times the search is (re)submitted if it fails on the NCBI side

# translation of the 'qblast'-style parameter names to the NCBI URL API ones
BLAST_PARAMETERS = {
    'composition_based_statistics': 'COMPOSITION_BASED_STATISTICS',
    'entrez_query': 'ENTREZ_QUERY',
    'expect': 'EXPECT',
    'filter': 'FILTER',
    'gapcosts': 'GAPCOSTS',
    'matrix_name': 'MATRIX_NAME',
    'megablast': 'MEGABLAST',
    'nucl_penalty': 'NUCL_PENALTY',
    'nucl_reward': 'NUCL_REWARD',
    'perc_ident': 'PERC_IDENT',
    'threshold': 'THRESHOLD',
    'word_size': 'WORD_SIZE',
}

def blast_extra_parameters(email:str,
                           tool:str,
                           **kwargs) -> dict:
    '''
    Build the extra parameters for the NCBI URL API from the given keyword arguments

    Parameters
    ----------
    email : str
        Contact email sent to NCBI
    tool : str
        Name of the tool sending the requests to NCBI
    **kwargs
        Additional parameters, either with the 'qblast'-style names (e.g. 'expect')
        or already with the NCBI URL API names (e.g. 'EXPECT')

    Returns
    -------
    dict
        Parameters ready to be sent to NCBI (excluding None values)
    '''
    parameters = {}
    for key, value in kwargs.items():
        if key.lower() in BLAST_PARAMETERS:
            key_url = BLAST_PARAMETERS[key.lower()]
        elif key.upper() in BLAST_PARAMETERS.values():
            key_url = key.upper()
        else:
            key_url = key.upper()
            logger.warning(f"WARNING: unknown BLAST parameter '{key}' - sending it to NCBI as '{key_url}'")
        parameters[key_url] = value
    parameters['EMAIL'] = email
    parameters['TOOL'] = tool
    return {key: value for key, value in parameters.items() if value is not None}

def blast_submit(seq_record:SeqRecord,
                 database:str,
                 hitlist:int,
                 parameters:dict) -> tuple:
    '''
    Send a BLAST search to NCBI with the URL API (CMD=Put)

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    database : str
        Database to search in
    hitlist : int
        Number of hits to request
    parameters : dict
        Extra parameters for the NCBI URL API

    Returns
    -------
    rid : str
        Request identifier (RID) of the submitted BLAST search
    rtoe : int
        Estimated wait time given by NCBI, in seconds (None if not provided)
    '''
    query_parameters = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': database,
        'QUERY': str(seq_record.seq),
        'HITLIST_SIZE': hitlist,
        'FORMAT_TYPE': 'XML',
    }
    query_parameters.update(parameters)
    logger.debug(f"BLAST query parameters: {query_parameters}")
    response = requests.post(NCBI_BLAST_URL, data=query_parameters, timeout=120)
    response.raise_for_status()
    rid = re.search(r'RID\s*=\s*(\S+)', response.text)
    rtoe = re.search(r'RTOE\s*=\s*(\d+)', response.text)
    if rid is None:
        logger.error(f"ERROR: no RID found in the NCBI response: {response.text[:500]}")
        raise RuntimeError("NCBI did not return a RID for the BLAST search")
    return rid.group(1), (int(rtoe.group(1)) if rtoe else None)

def blast_status(rid:str,
                 parameters:dict) -> tuple:
    '''
    Check the status of a BLAST search with the NCBI URL API (CMD=Get&FORMAT_OBJECT=SearchInfo)

    Parameters
    ----------
    rid : str
        Request identifier (RID) of the BLAST search
    parameters : dict
        Extra parameters for the NCBI URL API

    Returns
    -------
    status : str
        'WAITING', 'READY', 'FAILED' or 'UNKNOWN'
    has_hits : bool
        Whether the search found hits (None while the search is not ready)
    '''
    query_parameters = {'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid}
    query_parameters.update(parameters)
    response = requests.get(NCBI_BLAST_URL, params=query_parameters, timeout=60)
    response.raise_for_status()
    status = re.search(r'Status\s*=\s*(\w+)', response.text)
    if status is None:
        # an empty page means that the search is still running
        return 'WAITING', None
    has_hits = re.search(r'ThereAreHits\s*=\s*(\w+)', response.text)
    return status.group(1).upper(), (has_hits.group(1).lower() == 'yes' if has_hits else None)

def blast_wait(rid:str,
               parameters:dict,
               timeout:int,
               poll_interval:int,
               rtoe:int = None) -> str:
    '''
    Wait until a BLAST search is ready, reporting its status in the log

    Contrary to 'NCBIWWW.qblast', the status of the search is written in the log
    (so a long search can be told apart from a stuck one) and the search is not
    polled forever: after 'timeout' seconds the status 'TIMEOUT' is returned.

    Parameters
    ----------
    rid : str
        Request identifier (RID) of the BLAST search
    parameters : dict
        Extra parameters for the NCBI URL API
    timeout : int
        Maximum number of seconds to wait
    poll_interval : int
        Seconds between status checks
    rtoe : int, optional
        Estimated wait time given by NCBI, in seconds

    Returns
    -------
    str
        'READY', 'FAILED', 'UNKNOWN' or 'TIMEOUT'
    '''
    start_time = time.time()
    if rtoe:
        logger.debug(f"Waiting {min(rtoe, poll_interval)}s (NCBI estimated wait: {rtoe}s) before the first status check")
        time.sleep(min(rtoe, poll_interval))
    while True:
        elapsed = time.time() - start_time
        if elapsed >= timeout:
            logger.error(f"ERROR: BLAST search {rid} is taking too long (timeout: {timeout}s)")
            logger.error(f"ERROR: the results can still be retrieved later from: {NCBI_BLAST_URL}?CMD=Get&RID={rid}")
            return 'TIMEOUT'
        try:
            status, has_hits = blast_status(rid, parameters)
        except requests.RequestException as e:
            logger.warning(f"WARNING: could not check the status of BLAST {rid} ({e}). Retrying")
            status, has_hits = 'WAITING', None
        logger.info(f"-- BLAST {rid}: {status} (elapsed {elapsed/60:.1f} min){' - still searching...' if status == 'WAITING' else ''}")
        if status != 'WAITING':
            if status == 'READY' and has_hits is False:
                logger.warning("WARNING: the BLAST search finished with no hits")
            return status
        time.sleep(poll_interval)

def blast_download(rid:str,
                   parameters:dict) -> str:
    '''
    Download the results of a ready BLAST search with the NCBI URL API (CMD=Get)

    Parameters
    ----------
    rid : str
        Request identifier (RID) of the BLAST search
    parameters : dict
        Extra parameters for the NCBI URL API

    Returns
    -------
    str
        BLAST results in XML format
    '''
    query_parameters = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'}
    query_parameters.update(parameters)
    response = requests.get(NCBI_BLAST_URL, params=query_parameters, timeout=600)
    response.raise_for_status()
    if 'QBlastInfoBegin' in response.text or not response.text.lstrip().startswith('<?xml'):
        # while the search is not ready (or on error) the NCBI answers with a status page, not with XML results
        logger.error(f"ERROR: could not download the BLAST results for {rid}: {response.text[:500]}")
        raise RuntimeError(f"Could not download the BLAST results for {rid}")
    return response.text

def blast_remote(seq_record:SeqRecord,
                 blast_filename:str,
                 database:str,
                 hitlist:int,
                 rid:str = None,
                 email:str = NCBI_EMAIL,
                 tool:str = NCBI_TOOL,
                 timeout:int = BLAST_TIMEOUT,
                 poll_interval:int = BLAST_POLL_INTERVAL,
                 max_attempts:int = BLAST_MAX_ATTEMPTS,
                 **kwargs) -> None:
    '''
    Perform a remote BLAST search of a sequence with the NCBI URL API

    Unlike 'NCBIWWW.qblast', the status of the search is written in the log and
    searches failing on the NCBI side (status 'FAILED'/'UNKNOWN') are re-issued
    instead of being polled forever. The RID of the search is stored next to the
    results (so the results can be retrieved later) and a previous search can be
    resumed by passing its 'rid'.

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    database : str
        Database to search in
    hitlist : int
        Number of hits to request
    rid : str, optional
        RID of a previously submitted search to resume (def: submit a new search)
    email : str, optional
        Contact email sent to NCBI (def: 'NCBI_EMAIL')
    tool : str, optional
        Name of the tool sending the requests to NCBI (def: 'NCBI_TOOL')
    timeout : int, optional
        Maximum number of seconds to wait for each search (def: 'BLAST_TIMEOUT')
    poll_interval : int, optional
        Seconds between status checks (def: 'BLAST_POLL_INTERVAL')
    max_attempts : int, optional
        Times the search is (re)submitted if it fails on the NCBI side (def: 'BLAST_MAX_ATTEMPTS')
    **kwargs
        Additional parameters for the NCBI URL API (e.g. 'expect' or 'EXPECT')
    '''
    rid_filename = os.path.splitext(blast_filename)[0] + '.rid'
    parameters = blast_extra_parameters(email, tool, **kwargs)
    for attempt in range(1, max_attempts + 1):
        if rid is None:
            logger.info(f"-- Submitting BLAST search to NCBI (attempt {attempt}/{max_attempts})")
            rid, rtoe = blast_submit(seq_record, database, hitlist, parameters)
            with open(rid_filename, 'w') as f:
                f.write(rid)
            logger.info(f"-- BLAST request sent. RID: {rid} (NCBI estimated wait: {rtoe}s)")
            logger.info(f"-- Follow the search at: {NCBI_BLAST_URL}?CMD=Get&RID={rid}")
        else:
            logger.info(f"-- Resuming BLAST search. RID: {rid}")
            rtoe = None
        status = blast_wait(rid, parameters, timeout, poll_interval, rtoe)
        if status == 'READY':
            logger.info(f"-- BLAST {rid}: READY. Downloading the results")
            xml_results = blast_download(rid, parameters)
            with open(blast_filename, 'w') as f:
                f.write(xml_results)
            logger.info(f"-- BLAST results saved to '{blast_filename}' ({len(xml_results)/1e6:.1f} MB)")
            return
        if status == 'TIMEOUT':
            raise TimeoutError(f"BLAST search {rid} timed out after {timeout}s (its results can still be retrieved with RID '{rid}')")
        logger.error(f"ERROR: BLAST search {rid} ended with status '{status}'")
        if attempt == max_attempts:
            raise RuntimeError(f"BLAST search {rid} failed with status '{status}' after {max_attempts} attempts")
        logger.info(f"-- Re-issuing the BLAST search")
        rid = None
    # unreachable, but keeps linters happy
    raise RuntimeError("BLAST search failed")

def blastp_local(seq_record:SeqRecord,
                 blast_filename:str,
                 threads:int=0,
                 **kwargs) -> None:
    '''
    Run locally 'blastp' for a single sequence record against a database

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    threads : int, optional
        Number of threads to use (def: all available)
    '''
    threads = threads or len(os.sched_getaffinity(0))
    with NamedTemporaryFile(mode='w') as f:
        f.write(seq_record.format('fasta'))
        f.flush()
        additional_args = []
        for i, j in kwargs.items():
            additional_args.append(f'-{i}')
            additional_args.append(f'{j}')
        blastp_cmd = ['blastp', '-query', f.name, '-out', blast_filename, '-outfmt', '5', '-num_threads', str(threads)] + additional_args
        logger.debug(f"Running 'blastp' for {seq_record.id}:  {' '.join(blastp_cmd)}")
        blastp_run = subprocess.run(blastp_cmd)
        if blastp_run.returncode != 0:
            logger.error(f"ERROR: local 'blastp' for {seq_record.id} failed")
            raise RuntimeError(f"Error running 'blastp' for {seq_record.id}")

def blast_from_xml(seq_record:SeqRecord,
                   xml_path:str,
                   blast_filename:str) -> None:
    '''
    Use an existing BLAST XML file (e.g. downloaded from the NCBI web) instead of running a search
    '''
    if not os.path.isfile(xml_path):
        raise FileNotFoundError(f"BLAST XML file not found: '{xml_path}'")
    with open(xml_path, "r") as f:
        xml_string = f.read().replace('CREATE_VIEW', '')
    record = NCBIXML.read(StringIO(xml_string))
    query_len = getattr(record, 'query_letters', None)
    if query_len and query_len != len(seq_record.seq):
        raise ValueError(f"The XML '{xml_path}' is for a query of {query_len} aa, but {seq_record.id} has {len(seq_record.seq)} aa")
    logger.info(f"-- Skipping BLAST search: using existing results '{xml_path}' ({len(record.alignments)} alignments)")
    shutil.copyfile(xml_path, blast_filename)

def blast(seq_record:SeqRecord,
          blast_filename:str,
          database:str = None,
          hitlist:int = None,
          local:bool = False,
          **kwargs) -> None:
    '''
    Perform a BLAST search of a sequence

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    database : str, optional
        Database to search in (def: taken from 'parameters_default')
    hitlist : int, optional
        Number of hits to request (def: taken from 'parameters_default')
    local : bool, optional
        Use local BLASTp (def: False)
        BLAST+ must be installed and configured (executable in PATH and databases in BLASTDB)
    **kwargs
        Additional keyword arguments to pass to the 'blast_remote' function or 'blastp' program
    '''
    hitlist = hitlist or parameters_default['hitlist_size']
    database = database or parameters_default['blast_database']
    logger.info(dedent(f"""\
        -- Searching in BLAST ({'local' if local else 'remote'})
        Hitlist size:    {hitlist}
        Database:        {database}
        Output filename: '{blast_filename}'
        """))
    if not local:
        blast_remote(seq_record, blast_filename, database, hitlist, **kwargs)
    else:
        kwargs['db'] = database
        kwargs['max_target_seqs'] = hitlist
        blastp_local(seq_record, blast_filename, **kwargs)

def blast_parse(blast_filename:str,
                found_sequences_fastafile:str,
                min_identity:float,
                include_low_identity_seqs:bool = True,
                max_low_identity_seqs:int = -1,
                max_found_sequences:int = -1,
                lengths:list[int]=[]) -> int:
    '''
    Parse BLAST results to obtain the unique found sequences

    Read a BLAST's XML file with the results of a search, retrieve the
    HSP (high-scoring segment pairs) with a higher identity and write
    the whole corresponding sequences based on ther accession number.

    Parameters
    ----------
    blast_filename : str
        Path of the file containing BLAST results (XML format)
    found_sequences_fastafile : str
        Path of the file to write the found sequences (FASTA format)
    min_identity : float
        Minimum identity fraction expected in the found hits (0 to 1)
    include_low_identity_seqs : bool, optional
        Include low identity hits in the found sequences (def: True)
    max_low_identity_seqs : int, optional
        Maximum number of low identity hits, infinite by default (def: -1)
    max_found_sequences : int, optional
        Maximum number of found sequences, infinite by default (def: -1)
    lengths : list[int], optional
        Lengths of the sequences to restrict to, minimum and maximum (def: all)

    Returns
    -------
    n_low_identity_hsp : int
        Number of found hits with identity lower than the expected
    '''
    max_low_identity_seqs = max_low_identity_seqs if max_low_identity_seqs >= 0 else float('inf')
    max_found_sequences = max_found_sequences if max_found_sequences >= 0 else float('inf')
    logger.debug(f"Parsing BLAST results from '{blast_filename}'")
    with open(blast_filename, "r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    n_hsp = 0
    n_low_identity_hsp = 0
    accession_numbers = []
    if not records.alignments:
        logger.error(f"No BLAST alignments found in '{blast_filename}'")
        raise RuntimeError("No BLAST alignments found")
    for alignment in records.alignments:
        if 'partial' in alignment.title.lower():
            logger.debug(f"Partial found. Skipping sequence '{alignment.title}'")
            continue
        accession_number = alignment.accession
        logger.debug(f"\nAccession number: {accession_number}")
        for hsp in alignment.hsps:
            logger.debug(f"> {alignment.title}\n{hsp.sbjct}")
            identity = hsp.identities / alignment.length
            if identity < min_identity:
                logger.debug(f"Low identity hit found: {identity:.2f} ({'seq included' if include_low_identity_seqs else 'seq not included'})")
                if not include_low_identity_seqs:
                    continue
                n_low_identity_hsp += 1
            if include_low_identity_seqs or identity >= min_identity:
                accession_numbers.append(accession_number)
                n_hsp += 1
        if n_low_identity_hsp >= max_low_identity_seqs:
            include_low_identity_seqs = False
        if len(accession_numbers) >= max_found_sequences:
            logger.info("Maximum number of required found sequences reached")
            break
    accession_numbers = list(set(accession_numbers))
    logger.info(f"\nFound {len(accession_numbers)} sequences from {n_hsp} high-scoring segments ({n_low_identity_hsp} with low identity)\n")
    logger.debug(f"Fetching sequences with known accession numbers")
    get_seq(accession_numbers, found_sequences_fastafile, lengths=lengths, spaces=False)
    return n_low_identity_hsp

def blast_pro(seq_record:SeqRecord,
              blast_filename:str,
              found_sequences_fastafile:str,
              min_identity:float,
              min_low_identity_seqs:int,
              max_low_identity_seqs:int,
              additional_hits:int,
              hitlist:int,
              lengths:list[int]=[],
              include_low_identity_seqs:bool = True,
              database:str = None,
              local:bool = False,
              **kwargs) -> None:
    '''
    Perform a Pro BLAST search of a sequence

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    found_sequences_fastafile : str
        Path of the file to write the found sequences (FASTA format)
    min_identity : float
        Minimum identity fraction expected in the found identity (0 to 1)
    min_low_identity_seqs : int
        Minimum number of low identity hits
    max_low_identity_seqs : int
        Maximum number of low identity hits
    additional_hits : int
        Number of additional hits to add on each iteration
    hitlist : int
        Initial number of hits to search
    lengths : list[int], optional
        Lengths of the sequences to restrict to, minimum and maximum (def: all)
    include_low_identity_seqs : bool, optional
        Include low identity hits in the found sequences (def: True)
    database : str
        Database to search in (def: taken from 'parameters_default')
    local : bool, optional
        Use local BLASTp (def: False)
        BLAST+ must be installed and configured (executable in PATH and databases in BLASTDB)
    **kwargs
        Additional keyword arguments to pass to the 'blast_remote' function or 'blastp' program
    '''
    max_iter = 100
    for iteration in range(max_iter):
        logger.info(f"Pro BLAST iteration {iteration + 1}\n")
        blast(seq_record, blast_filename, database, hitlist, local, **kwargs)
        n_low_identity_seqs = blast_parse(blast_filename, found_sequences_fastafile, min_identity, include_low_identity_seqs, max_low_identity_seqs, -1, lengths)
        if n_low_identity_seqs >= min_low_identity_seqs:
            break
        hitlist += additional_hits
        logger.info(f"The number of low identity sequences is below the desired value\nBlasting again with {hitlist} hits")
        logger.debug(f"Removing BLAST results file '{blast_filename}' and found sequences file '{found_sequences_fastafile}'")
        os.remove(blast_filename)
        os.remove(found_sequences_fastafile)
    else:
        logger.error("ERROR: Pro BLAST failed: maximum number of iterations reached")
        raise Exception("Maximum number of iterations reached")
    logger.info(f"\nPro BLAST done succesfully\n")
