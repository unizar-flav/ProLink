
import logging

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def get_seq(IDs:list[str], fastafile:str='', lengths:list[int]=[], spaces:bool=True) -> list[SeqRecord]:
    '''
    Obtain the sequences of protein's IDs

    Parameters
    ----------
    IDs : list[str]
        List of protein's IDs
    fastafile : str, optional
        Path of a fasta file to write the sequences
    lengths : list[int], optional
        Lengths of the sequences to restrict to, minimum and maximum (def: all)
    spaces : bool, optional
        Include spaces in the descriptions of the FASTA file
        Otherwise, replace them with underscores (def: True)

    Returns
    -------
    list[SeqRecord]
        List of SeqRecord objects corresponding to the sequences of the IDs
    '''
    seq_list = []
    IDs = [IDs] if isinstance(IDs, str) else IDs
    Entrez.email = "example@example.com"
    with Entrez.efetch(db="protein", rettype="gb", retmode="text", id=IDs) as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            rec = SeqRecord(Seq(str(seq_record.seq)), id=seq_record.id, description=seq_record.description)
            seq_list.append(rec)
    if len(seq_list) != len(IDs):
        not_found = set(IDs) - set([seq_record.id for seq_record in seq_list])
        logger.warning(f"WARNING: Unable to obtain sequences for: {', '.join(not_found)}")
    if lengths:
        logger.debug(f"Discarting sequences not in the range {lengths}")
        for seq_record in seq_list:
            if not (lengths[0] <= len(seq_record.seq) <= lengths[1]):
                logger.warning(f"WARNING: Discarting sequence '{seq_record.id}' with length {len(seq_record.seq)} {'+' if lengths[0] < len(seq_record.seq) else '-'}")
        seq_list = [seq_record for seq_record in seq_list if lengths[0] <= len(seq_record.seq) <= lengths[1]]
    if not spaces:
        logger.debug("Replacing spaces with underscores in the descriptions of the FASTA file")
        for seq_record in seq_list:
            seq_record.description = seq_record.description.replace(" ", "_")
    if fastafile:
        logger.debug(f"Saving sequences to '{fastafile}'")
        SeqIO.write(seq_list, fastafile, "fasta")
    return seq_list

def check_seq_in(seq:SeqRecord, fastafile:str, rewrite:bool=True, spaces:bool=True) -> bool:
    '''
    Check if a sequence is in a FASTA file

    Parameters
    ----------
    seq : SeqRecord
        Sequence to check
    fastafile : str
        Path of the fasta file to read (and write if not found)
    rewrite : bool, optional
        Rewrite the fasta file with the query sequence if not found (def: True)
    spaces : bool, optional
        Include spaces in the description of the sequence in the FASTA file
        Otherwise, replace them with underscores (def: True)

    Returns
    -------
    bool
        sequence in the FASTA file or not
    '''
    fastafile_seqs = list(SeqIO.parse(fastafile, "fasta"))
    seq_in = any([seq.seq == fastafile_seq.seq for fastafile_seq in fastafile_seqs])
    logger.info(f"Query sequence {'found' if seq_in else 'not found'} in '{fastafile}'")
    if not seq_in and rewrite:
        logger.info(f"Prepending query sequence in '{fastafile}'")
        if not spaces:
            logger.debug("Replacing spaces with underscores in the description")
            seq.description = seq.description.replace(" ", "_")
        fastafile_seqs.insert(0, seq)
        SeqIO.write(fastafile_seqs, fastafile, "fasta")
    return seq_in
