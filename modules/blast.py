from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from io import StringIO
import subprocess
import os

def blast(hitlist_range, blast_database, blast_filename,  my_seq_record):
    print("----Searching in BLAST----")
    print("Hitlist size:    " + str(hitlist_range))
    print("Database:        " + blast_database)
    print("Output filename: " + blast_filename)
    
    result_handle = NCBIWWW.qblast("blastp", str(blast_database), my_seq_record.seq, hitlist_size=hitlist_range)
    print(result_handle)

    with open(blast_filename, 'w') as save_file: 
        blast_results = result_handle.read()
        save_file.write(blast_results)

def parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, outputs_dir):
    with open(blast_filename,"r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
        low_identity_seqs = 0
    found_sequences=[]
    for alignment in records.alignments:
        for hsp in alignment.hsps:
            rec_f = SeqRecord(
                        Seq(hsp.sbjct,
                        ),
                        id=alignment.title,
                    )
            print('>', alignment.title) 
            print(hsp.sbjct)
            if (hsp.identities / alignment.length) < expected_min_identity:
                    low_identity_seqs += 1
                    print("Low identity sequence")
            print()
            rec_f.description = ""
            rec_f.id=rec_f.id.replace(" ", "_")
            rec_f.seq=rec_f.seq.ungap("-")
            found_sequences.append(rec_f)
    found_sequences_fastafile= "./" + outputs_dir + "/protein_" + str(my_sequence_index) + "/found_sequences.fasta"
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    print(low_identity_seqs)
    return low_identity_seqs

def s_parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, low_identity_seqs, max_low_identity_seqs, outputs_dir):
    with open(blast_filename,"r") as f:
        xml_string = f.read()
        xml_string = xml_string.replace('CREATE_VIEW', '')
        records = NCBIXML.read(StringIO(xml_string))
    found_sequences=[]
    sequence_index= 0
    for alignment in records.alignments:
        hsp = alignment.hsps[0]
        if sequence_index < 10000:
          if low_identity_seqs < max_low_identity_seqs:
              rec_f = SeqRecord(
                          Seq(hsp.sbjct,
                          ),
                          id=alignment.title,
                      )
              sequence_index += 1 
              print("Sequence num " + str(sequence_index))
              print('>', alignment.title) 
              print(hsp.sbjct)
              print(expected_min_identity)
              if (hsp.identities / alignment.length) < expected_min_identity:
                  low_identity_seqs += 1
                  print("low identity seq!")
              print()
              rec_f.description = ""
              rec_f.id=rec_f.id.replace(" ", "_")
              rec_f.seq=rec_f.seq.ungap("-")
              found_sequences.append(rec_f)
          else:
              break
        else:
          break
    found_sequences_fastafile= "./" + outputs_dir + "/protein_" + str(my_sequence_index) + "/found_sequences.fasta"
    SeqIO.write(found_sequences, found_sequences_fastafile, "fasta")
    print(sequence_index)
    print(low_identity_seqs)
    return low_identity_seqs
  
def s_blast(my_sequence_index, blast_database, hitlist_range, my_seq_record, blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, outputs_dir):  
    blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
    
    if parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, outputs_dir) < min_low_identity_seqs:
        print()
        print("The number of low identity sequences is below the desired value")
        os.remove(blast_filename)
        os.remove(found_sequences_fastafile)
        hitlist_range += additional_hits
        print("Blasting again with " + str(hitlist_range) + " hits")
        blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
        low_identity_seqs = 0
        while s_parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, low_identity_seqs, max_low_identity_seqs, outputs_dir) < min_low_identity_seqs:
            print()
            print("The number of low identity sequences is below the desired value")
            os.remove(blast_filename)
            os.remove(found_sequences_fastafile)
            hitlist_range += additional_hits
            print("Blasting again with " + str(hitlist_range) + " hits")
            blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
            low_identity_seqs = 0
            s_parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, low_identity_seqs, max_low_identity_seqs, outputs_dir)
    print()
    print("Smart blast done succesfully")
    print()
