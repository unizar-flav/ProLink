
import logging

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def get_seq(IDs:list[str], fastafile:str=None) -> list[SeqRecord]:
    '''
    Obtain the sequences of protein's IDs

    Parameters
    ----------
    IDs : list[str]
        List of protein's IDs
    fastafile : str, optional
        Path of a fasta file to write the sequences

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
    if fastafile:
        logger.debug(f"Saving sequences to '{fastafile}'")
        SeqIO.write(seq_list, fastafile, "fasta")
    return seq_list
