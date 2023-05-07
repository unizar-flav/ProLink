
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    print("Obtaining sequences\n")
    seq_list = []
    Entrez.email = "example@example.com"
    with Entrez.efetch(db="protein", rettype="gb", retmode="text", id=IDs) as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            rec = SeqRecord(Seq(str(seq_record.seq)), id=seq_record.id, description=seq_record.description)
            seq_list.append(rec)
            print(f"> {seq_record.id} {seq_record.description}\n{seq_record.seq}\n")
    if fastafile:
        SeqIO.write(seq_list, fastafile, "fasta")
    return seq_list
