
import os
from io import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .subprocess_functions import blastp
from .. import parameters_default


def blast(seq_record:SeqRecord, blast_filename:str, database:str=None, hitlist:int=None, local:bool=False, **kwargs) -> None:
    '''
    Perform a BLAST search of a sequence

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    database : str, optional
        Database to search in (def: taken from 'parameters.cfg')
    hitlist : int, optional
        Number of hits to request (def: taken from 'parameters.cfg')
    local : bool, optional
        Use local BLASTp (def: False)
        BLAST+ must be installed and configured (executable in PATH and databases in BLASTDB)
    **kwargs
        Additional keyword arguments to pass to the 'qblast' function or 'blastp' program
    '''
    if not local:
        kwargs['program'] = 'blastp'
        kwargs['database'] = database or parameters_default['blast_database']
        kwargs['sequence'] = seq_record.seq
        kwargs['hitlist_size'] = hitlist or parameters_default['hitlist_size']
        kwargs['format_type'] = 'XML'
        print(f"---- Searching in BLAST (remote) ----")
        print(f"Hitlist size:    {kwargs['hitlist_size']}")
        print(f"Database:        {kwargs['database']}")
        print(f"Output filename: {blast_filename}")
        result_handle = NCBIWWW.qblast(**kwargs)
        with open(blast_filename, 'w') as f:
            f.write(result_handle.read())
    else:
        kwargs['db'] = database or parameters_default['blast_database']
        kwargs['max_target_seqs'] = hitlist or parameters_default['hitlist_size']
        blastp(seq_record, blast_filename, **kwargs)


def blast_parse(blast_filename:str, found_sequences_fastafile:str, expected_min_identity:float, remove_gaps:bool=True, include_low_identity:bool=True, max_low_identity_seqs:int=None, max_found_sequences:int=None) -> int:
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
    with open(blast_filename,"r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    n_low_identity_seqs = 0
    found_sequences=[]
    for alignment in records.alignments:
        for hsp in alignment.hsps:
            rec_f = SeqRecord(Seq(hsp.sbjct), id=alignment.title)
            print('>', alignment.title)
            print(hsp.sbjct)
            identity = hsp.identities / alignment.length
            if identity < expected_min_identity:
                n_low_identity_seqs += 1
                print("Low identity sequence!")
            print()
            rec_f.description = ""
            rec_f.id=rec_f.id.replace(" ", "_").replace(" <unknown description>", "")
            if remove_gaps:
                rec_f.seq=rec_f.seq.replace("-", "")
            if include_low_identity or identity >= expected_min_identity:
                found_sequences.append(rec_f)
        if max_found_sequences and len(found_sequences) >= max_found_sequences or \
           max_low_identity_seqs and n_low_identity_seqs >= max_low_identity_seqs:
            break
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    if not include_low_identity:
        n_low_identity_seqs = 0
    return n_low_identity_seqs


def p_blast(seq_record:SeqRecord, blast_filename:str, found_sequences_fastafile:str, remove_gaps:bool, expected_min_identity:int, min_low_identity_seqs:int, max_low_identity_seqs:int, additional_hits:int, hitlist:int, database:str=None, local:bool=False, **kwargs) -> None:
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
    expected_min_identity : int
        Minimum identity percentage expected in the found sequences
    min_low_identity_seqs : int
        Minimum number of low identity sequences
    max_low_identity_seqs : int
        Maximum number of low identity sequences
    additional_hits : int
        Number of additional hits to add on each iteration
    hitlist : int
        Initial number of hits to search
    database : str
        Database to search in (def: taken from 'parameters.cfg')
    local : bool, optional
        Use local BLASTp (def: False)
        BLAST+ must be installed and configured (executable in PATH and databases in BLASTDB)
    **kwargs
        Additional keyword arguments to pass to the 'qblast' function or 'blastp' program
    '''
    include_low_identity = True
    blast(seq_record, blast_filename, database, hitlist, local, **kwargs)
    if blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, remove_gaps, include_low_identity) < min_low_identity_seqs:
        print("\nThe number of low identity sequences is below the desired value")
        os.remove(blast_filename)
        os.remove(found_sequences_fastafile)
        hitlist += additional_hits
        print(f"Blasting again with {hitlist} hits")
        blast(seq_record, blast_filename, database, hitlist, local, **kwargs)
        while blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, remove_gaps, include_low_identity, max_low_identity_seqs, 10000) < min_low_identity_seqs and hitlist < 10000:
            print("\nThe number of low identity sequences is below the desired value")
            os.remove(blast_filename)
            os.remove(found_sequences_fastafile)
            hitlist += additional_hits
            blast(seq_record, blast_filename, database, hitlist, local, **kwargs)
            blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, remove_gaps, include_low_identity, max_low_identity_seqs, 10000)
    print("\nPro blast done succesfully\n")
