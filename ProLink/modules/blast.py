
import logging
import os
from io import StringIO
from textwrap import dedent

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .. import parameters_default
from .subprocess_functions import blastp_run


logger = logging.getLogger()

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
        blastp_run(seq_record, blast_filename, **kwargs)

def blast_parse(blast_filename:str,
                found_sequences_fastafile:str,
                expected_min_identity:float,
                remove_gaps:bool = True,
                include_low_identity:bool = True,
                max_low_identity_seqs:int = None,
                max_found_sequences:int = None) -> int:
    '''
    Parse BLAST results

    Parameters
    ----------
    blast_filename : str
        Path of the file containing BLAST results (XML format)
    found_sequences_fastafile : str
        Path of the file to write the found sequences (FASTA format)
    expected_min_identity : float
        Minimum identity percentage expected in the found sequences (0 to 1)
    remove_gaps : bool, optional
        Remove gaps from the found sequences (def: True)
    include_low_identity : bool, optional
        Include low identity sequences in the found sequences (def: True)
    max_low_identity_seqs : int, optional
        Maximum number of low identity sequences, infinite by default (def: None)
    max_found_sequences : int, optional
        Maximum number of found sequences, infinite by default (def: None)

    Returns
    -------
    n_low_identity_seqs : int
        Number of found sequences with identity lower than the expected
    '''
    logger.debug(f"Parsing BLAST results from '{blast_filename}'")
    with open(blast_filename, "r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    n_low_identity_seqs = 0
    found_sequences = []
    for alignment in records.alignments:
        if 'partial' in alignment.title.lower():
            logger.debug(f"Skipping partial sequence '{alignment.title}'")
            continue
        for hsp in alignment.hsps:
            rec_f = SeqRecord(Seq(hsp.sbjct), id=alignment.title)
            logger.debug(f"\n> {alignment.title}\n{hsp.sbjct}")
            identity = hsp.identities / alignment.length
            if identity < expected_min_identity:
                logger.info(f"Low identity sequence found: {identity:.2f} ({'included' if include_low_identity else 'not included'})")
                if not include_low_identity:
                    continue
                n_low_identity_seqs += 1
            rec_f.description = ""
            rec_f.id = rec_f.id.replace(" <unknown description>", "").replace(" ", "_")
            if remove_gaps:
                rec_f.seq = rec_f.seq.replace("-", "")
            if include_low_identity or identity >= expected_min_identity:
                found_sequences.append(rec_f)
        if max_low_identity_seqs and n_low_identity_seqs >= max_low_identity_seqs:
            include_low_identity = False
        if max_found_sequences and len(found_sequences) >= max_found_sequences:
            logger.info("Maximum number of required found sequences reached")
            break
    logger.debug(f"Writing found sequences to '{found_sequences_fastafile}'")
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    logger.info(f"Found {len(found_sequences)} sequences ({n_low_identity_seqs} with low identity)\n")
    return n_low_identity_seqs

def p_blast(seq_record:SeqRecord,
            blast_filename:str,
            found_sequences_fastafile:str,
            remove_gaps:bool,
            expected_min_identity:float,
            min_low_identity_seqs:int,
            max_low_identity_seqs:int,
            additional_hits:int,
            hitlist:int,
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
    remove_gaps : bool
        Remove gaps from the found sequences
    expected_min_identity : float
        Minimum identity percentage expected in the found sequences (0 to 1)
    min_low_identity_seqs : int
        Minimum number of low identity sequences
    max_low_identity_seqs : int
        Maximum number of low identity sequences
    additional_hits : int
        Number of additional hits to add on each iteration
    hitlist : int
        Initial number of hits to search
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
        n_low_identity_seqs = blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, remove_gaps, True, max_low_identity_seqs, max_hitlist)
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
