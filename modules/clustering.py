from pprint import pformat
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from csv import reader
from ProLink.modules.pfam import search_hmmer_pfam
import os

def cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile):
    print("Clustering sequences with ALFATClust")
    #subprocess.call(['alfatclust.py', '-i', found_sequences_fastafile, '-o', cluster_results_file, '-e', cluster_evaluation_file, '-ea' ,'-l' , str(similarity)])
    cmd_line = 'alfatclust.py'+' -i '+ str(found_sequences_fastafile)+' -o '+str(cluster_results_file)+' -e '+str(cluster_evaluation_file)+' -ea '+'-l '+str(similarity)
    ! {cmd_line}
    
    clustered_sequences = []
    number_of_clusters = 0
    my_sequence_domains = search_hmmer_pfam(str(my_seq_record.seq)).keys()
    read_obj = open(cluster_evaluation_file, 'r')
    csv_reader = reader(read_obj)
    header = next(csv_reader)
    if header != None:
        for row in csv_reader:
            cluster_center_seq_description = row[4]
            cluster_id = "Cluster:"+row[0] +"_Seqs:"+ row[1] +"_ASI:"+ row[2] +"_MSI:"+ row[3]
            for seq_record in SeqIO.parse(found_sequences_fastafile, "fasta"):
                if seq_record.description == cluster_center_seq_description:
                    rec_c = SeqRecord(
                        Seq(str(seq_record.seq)),
                        id= cluster_id,
                        description = cluster_center_seq_description.replace(" <unknown description>", "")
                    )
                    try:
                        subject_sequence_domains = search_hmmer_pfam(str(seq_record.seq)).keys()
                        print(subject_sequence_domains)
                    except KeyError:
                        print("No domains found")
                        subject_sequence_domains = "No domains found"
                    if my_sequence_domains != subject_sequence_domains:
                        domains = " Different_domains:" + str(subject_sequence_domains).replace("dict_keys", "").replace("([","").replace("])","").replace("'","")
                    else:
                        domains = " Same_domains:" + str(subject_sequence_domains).replace("dict_keys", "").replace("([","").replace("])","").replace("'","")
                    rec_c_m = SeqRecord(
                        Seq(str(seq_record.seq)),
                        id= cluster_id,
                        description = cluster_center_seq_description.replace(" <unknown description>", "") + domains
                    )
                    clustered_sequences.append(rec_c_m)
                    number_of_clusters += 1
                    print("Cluster number " + str(number_of_clusters))
                    print(seq_record.description)
                    print(seq_record.seq)
                    print()  
                    break
    SeqIO.write(clustered_sequences, cluster_results_fastafile, "fasta")
    print("Number of clusters: "+ str(number_of_clusters))
    return number_of_clusters



def s_cluster(found_sequences_fastafile, my_seq_record, similarity, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile):
    num_of_clusters = cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
    #print("Number of clusters: " + str(number_of_clusters))
    #print(str(cluster_results_file))
    while num_of_clusters <= min_number_of_clusters_to_cluster_again or num_of_clusters >= max_number_of_clusters_to_cluster_again:
        
        if min_number_of_clusters_to_cluster_again < num_of_clusters:
            print("The number of clusters is between the desired values")
            print("Smart clustering done succesfully")

        if num_of_clusters < min_number_of_clusters_to_cluster_again:
            print("The number of clusters is below the minimum")
            os.remove(cluster_results_file)
            os.remove(cluster_evaluation_file)
            os.remove(cluster_results_fastafile)
            similarity += 0.1
            print("Clustering again with " + str(similarity) + " similarity")
            num_of_clusters = cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)

        elif num_of_clusters >= max_number_of_clusters_to_cluster_again:
            print("The number of clusters exceeds the maximum")
            os.remove(cluster_results_file)
            os.remove(cluster_evaluation_file)
            os.remove(cluster_results_fastafile)
            similarity -= 0.1
            print("Clustering again with " + str(similarity) + " similarity")
            num_of_clusters = cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
