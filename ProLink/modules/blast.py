
import logging
import os
import subprocess
from io import StringIO
from multiprocessing import cpu_count
from tempfile import NamedTemporaryFile
from textwrap import dedent

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord

from .. import parameters_default
from .obtaining_sequences import get_seq


logger = logging.getLogger()

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
    threads = threads or cpu_count()
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
        Additional keyword arguments to pass to the 'qblast' function or 'blastp' program
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
        kwargs['program'] = 'blastp'
        kwargs['database'] = database
        kwargs['sequence'] = seq_record.seq
        kwargs['hitlist_size'] = hitlist
        kwargs['format_type'] = 'XML'
        logger.debug(f"BLAST query parameters: {kwargs}")
        result_handle = NCBIWWW.qblast(**kwargs)
        with open(blast_filename, 'w') as f:
            f.write(result_handle.read())
    else:
        kwargs['db'] = database
        kwargs['max_target_seqs'] = hitlist
        blastp_local(seq_record, blast_filename, **kwargs)

def blast_parse(blast_filename:str,
                found_sequences_fastafile:str,
                expected_min_identity:float,
                include_low_identity:bool = True,
                max_low_identity_seqs:int = None,
                max_found_sequences:int = None,
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
    expected_min_identity : float
        Minimum identity percentage expected in the found hits (0 to 1)
    include_low_identity : bool, optional
        Include low identity hits in the found sequences (def: True)
    max_low_identity_seqs : int, optional
        Maximum number of low identity hits, infinite by default (def: None)
    max_found_sequences : int, optional
        Maximum number of found sequences, infinite by default (def: None)
    lengths : list[int], optional
        Lengths of the sequences to restrict to, minimum and maximum (def: all)

    Returns
    -------
    n_low_identity_hsp : int
        Number of found hits with identity lower than the expected
    '''
    logger.debug(f"Parsing BLAST results from '{blast_filename}'")
    with open(blast_filename, "r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    n_hsp = 0
    n_low_identity_hsp = 0
    accession_numbers = []
    for alignment in records.alignments:
        if 'partial' in alignment.title.lower():
            logger.debug(f"Partial found. Skipping sequence '{alignment.title}'")
            continue
        accession_number = alignment.accession
        logger.debug(f"\nAccession number: {accession_number}")
        for hsp in alignment.hsps:
            logger.debug(f"> {alignment.title}\n{hsp.sbjct}")
            identity = hsp.identities / alignment.length
            if identity < expected_min_identity:
                logger.debug(f"Low identity hit found: {identity:.2f} ({'seq included' if include_low_identity else 'seq not included'})")
                if not include_low_identity:
                    continue
                n_low_identity_hsp += 1
            if include_low_identity or identity >= expected_min_identity:
                accession_numbers.append(accession_number)
                n_hsp += 1
        if max_low_identity_seqs and n_low_identity_hsp >= max_low_identity_seqs:
            include_low_identity = False
        if max_found_sequences and len(accession_numbers) >= max_found_sequences:
            logger.info("Maximum number of required found sequences reached")
            break
    accession_numbers = list(set(accession_numbers))
    logger.info(f"\nFound {len(accession_numbers)} sequences from {n_hsp} high-scoring segments ({n_low_identity_hsp} with low identity)\n")
    logger.debug(f"Fetching sequences with known accession numbers")
    get_seq(accession_numbers, found_sequences_fastafile, lengths=lengths, spaces=False)
    return n_low_identity_hsp

def p_blast(seq_record:SeqRecord,
            blast_filename:str,
            found_sequences_fastafile:str,
            expected_min_identity:float,
            min_low_identity_seqs:int,
            max_low_identity_seqs:int,
            additional_hits:int,
            hitlist:int,
            lengths:list[int]=[],
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
    expected_min_identity : float
        Minimum identity percentage expected in the found identity (0 to 1)
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
    database : str
        Database to search in (def: taken from 'parameters_default')
    local : bool, optional
        Use local BLASTp (def: False)
        BLAST+ must be installed and configured (executable in PATH and databases in BLASTDB)
    **kwargs
        Additional keyword arguments to pass to the 'qblast' function or 'blastp' program
    '''
    max_hitlist = 10000     # Google Colab limit
    max_iter = 100
    for iteration in range(max_iter):
        logger.info(f"Pro BLAST iteration {iteration + 1}\n")
        blast(seq_record, blast_filename, database, hitlist, local, **kwargs)
        n_low_identity_seqs = blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, True, max_low_identity_seqs, max_hitlist, lengths)
        if n_low_identity_seqs >= min_low_identity_seqs or hitlist >= max_hitlist:
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
