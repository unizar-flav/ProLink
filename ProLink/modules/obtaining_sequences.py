
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def obtain_fasta_file(id, outputs_dir):
    print("Obtaining sequences")
    print()
    seq_list = []
    Entrez.email = "A.N.Other@example.com"
    with Entrez.efetch(db="protein", rettype="gb", retmode="text", id=id) as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            rec = SeqRecord(
                Seq(str(seq_record.seq),
                ),
                id=seq_record.id,
                description=seq_record.description,
             )
            seq_list.append(rec)
            print(">%s %s" % (seq_record.id, seq_record.description))
            print(seq_record.seq)
            print()
        SeqIO.write(seq_list, "./" + outputs_dir +"/my_sequences.fasta", "fasta")
    return seq_list
