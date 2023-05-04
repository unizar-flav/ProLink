
import os
from io import StringIO

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .. import parameters_default


def blast(seq_record:SeqRecord, blast_filename:str, **qblast_kwargs) -> None:
    '''
    Perform a BLAST search of a sequence

    Parameters
    ----------
    seq_record : SeqRecord
        Sequence to search
    blast_filename : str
        Path of the file to write the BLAST results (XML format)
    **qblast_kwargs
        Keyword arguments to pass to the 'qblast' function
        i.e. database, hitlist_size, etc.
    '''
    # default needed values
    qblast_kwargs['program'] = 'blastp'
    qblast_kwargs['database'] = qblast_kwargs.get('database', parameters_default['blast_database'])
    qblast_kwargs['sequence'] = seq_record.seq
    qblast_kwargs['hitlist_size'] = qblast_kwargs.get('hitlist_size', parameters_default['hitlist_range'])
    qblast_kwargs['format_type'] = 'XML'
    print(f"---- Searching in BLAST ----")
    print(f"Hitlist size:    {qblast_kwargs['hitlist_size']}")
    print(f"Database:        {qblast_kwargs['database']}")
    print(f"Output filename: {blast_filename}")
    result_handle = NCBIWWW.qblast(**qblast_kwargs)
    with open(blast_filename, 'w') as f:
        f.write(result_handle.read())


def parse(blast_filename:str, found_sequences_fastafile:str, remove_gaps:bool, expected_min_identity:int) -> int:
    '''
    Parse BLAST results

    Parameters
    ----------
    blast_filename : str
        Path of the file containing BLAST results (XML format)
    found_sequences_fastafile : str
        Path of the file to write the found sequences (FASTA format)
    remove_gaps : bool
        Remove gaps from the found sequences
    expected_min_identity : int
        Minimum identity percentage expected in the found sequences

    Returns
    -------
    int
        Number of found sequences with identity lower than the expected
    '''
    with open(blast_filename,"r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    low_identity_seqs = 0
    found_sequences=[]
    for alignment in records.alignments:
        for hsp in alignment.hsps:
            rec_f = SeqRecord(Seq(hsp.sbjct), id=alignment.title)
            print('>', alignment.title)
            print(hsp.sbjct)
            if (hsp.identities / alignment.length) < expected_min_identity:
                    low_identity_seqs += 1
                    print("Low identity sequence")
            print()
            rec_f.description = ""
            rec_f.id=rec_f.id.replace(" ", "_").replace(" <unknown description>", "")
            if remove_gaps:
                rec_f.seq=rec_f.seq.replace("-", "")
            found_sequences.append(rec_f)
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    return low_identity_seqs


def p_parse(blast_filename:str, found_sequences_fastafile:str, remove_gaps:bool, expected_min_identity:int, min_low_identity_seqs:int, max_low_identity_seqs:int) -> int:
    '''
    Parse Pro BLAST results

    Parameters
    ----------
    blast_filename : str
        Path of the file containing BLAST results (XML format)
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

    Returns
    -------
    int
        Number of found sequences with identity lower than the expected
    '''
    with open(blast_filename,"r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    low_identity_seqs = min_low_identity_seqs
    found_sequences=[]
    sequence_index= 0
    for alignment in records.alignments:
        hsp = alignment.hsps[0]
        if sequence_index < 10000:
          if low_identity_seqs < max_low_identity_seqs:
              rec_f = SeqRecord(Seq(hsp.sbjct), id=alignment.title)
              sequence_index += 1
              print("Sequence num " + str(sequence_index))
              print('>', alignment.title)
              print(hsp.sbjct)
              if (hsp.identities / alignment.length) < expected_min_identity:
                  low_identity_seqs += 1
                  print("low identity seq!")
              print()
              rec_f.description = ""
              rec_f.id=rec_f.id.replace(" ", "_").replace(" <unknown description>", "")
              if remove_gaps:
                rec_f.seq=rec_f.seq.replace("-", "")
              found_sequences.append(rec_f)
          else:
              break
        else:
          break
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    return low_identity_seqs


def p_blast(seq_record:SeqRecord, blast_filename:str, found_sequences_fastafile:str, remove_gaps:bool, expected_min_identity:int, min_low_identity_seqs:int, max_low_identity_seqs:int, additional_hits:int, hitlist_range:int, **qblast_kwargs) -> None:
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
    hitlist_range : int
        Initial number of hits to search
    **qblast_kwargs
        Keyword arguments for the 'qblast' function
        i.e. database, hitlist_size, etc.
    '''
    blast(seq_record, blast_filename, hitlist_size=hitlist_range, **qblast_kwargs)
    if parse(blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity) < min_low_identity_seqs:
        print("\nThe number of low identity sequences is below the desired value")
        os.remove(blast_filename)
        os.remove(found_sequences_fastafile)
        hitlist_range += additional_hits
        print(f"Blasting again with {hitlist_range} hits")
        blast(seq_record, blast_filename, hitlist_size=hitlist_range, **qblast_kwargs)
        low_identity_seqs = 0
        while p_parse(blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, low_identity_seqs, max_low_identity_seqs) < min_low_identity_seqs and hitlist_range < 10000:
            print("\nThe number of low identity sequences is below the desired value")
            os.remove(blast_filename)
            os.remove(found_sequences_fastafile)
            hitlist_range += additional_hits
            blast(seq_record, blast_filename, hitlist_size=hitlist_range, **qblast_kwargs)
            low_identity_seqs = 0
            p_parse(blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, low_identity_seqs, max_low_identity_seqs)
    print("\nPro blast done succesfully\n")
