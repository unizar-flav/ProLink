from genericpath import isdir
import os
import subprocess
from .modules.obtaining_sequences import obtain_fasta_file
from .modules.blast import blast
from .modules.blast import parse
from .modules.blast import s_blast
from .modules.clustering import cluster
from .modules.clustering import s_cluster
from .modules.subprocess_functions import align
from .modules.subprocess_functions import weblogo3
from .modules.subprocess_functions import tree
import os
import sys
if os.path.isdir("./outputs"):
    print("outputs directory already exists")
else:
    os.mkdir("outputs")

#--------------------PARAMETERS----------------------
#For local executing of the script:
"""
#BLAST

query_proteins = "ABQ62066.1, ABQ62091.1, ABQ62490.1" #Write the protein IDs here. Eg: "ABQ62066.1, ABQ62091.1, ABQ62490.1"
hitlist_range = 5000 #Desired minimum number of found sequences.
blast_database = "refseq_protein"
smart_blast_ = True
#If smart_blast:
max_low_identity_seqs = 1
min_low_identity_seqs = 1 #If the number of sequences with a low identity identity percentage is below this value, the script will blast again until this value is reached.
expected_min_identity = 0.25 #Maximum identity percentage to consider a sequence a "low identity sequence".
additional_hits = 2000 #The new blast will search for hitlist_range + aditional_hits sequences.

remove_gaps = True

#CLUSTERING
cluster_seqs = True
similarity = 0.5 #Initial similarity treshold to group the sequences into clusters.
smart_clustering_ = True
#If smart_clustering:
min_number_of_clusters_to_cluster_again = 250
max_number_of_clusters_to_cluster_again = 700

#ALIGNMENT
align_seqs = True

#SEQUENCE LOGO GENERATION
generate_logo = True
#if logo:
weblogo_format = "png"

#PHYLOGENETIC TREE GENERATION
generate_tree = True
#If tree:
tree_type = "NJ" #Either write NJ (Neighbor joining) or ML (Maximum likehood).
bootstrap_replications = "500" #Write 250, 500, 1000, 2000 or 5000.

#For Google Colab executing:
query_proteins = str(sys.argv[1])
hitlist_range = int(sys.argv[2])
blast_database = str(sys.argv[3])
smart_blast_ = sys.argv[4]
#If smart_blast:
max_low_identity_seqs = int(sys.argv[5])
min_low_identity_seqs = int(sys.argv[6])
expected_min_identity = float(sys.argv[7])
additional_hits = int(sys.argv[8])

remove_gaps = sys.argv[9]

#CLUSTERING
cluster_seqs = sys.argv[20]
similarity = float(sys.argv[10])
smart_clustering_ = sys.argv[11]
#If smart_clustering:
min_number_of_clusters_to_cluster_again = int(sys.argv[12])
max_number_of_clusters_to_cluster_again = int(sys.argv[13])

#ALIGNMENT
align_seqs = sys.argv[14]

#SEQUENCE LOGO GENERATION
generate_logo = sys.argv[15]
#if logo:
weblogo_format = str(sys.argv[16])

#PHYLOGENETIC TREE GENERATION
generate_tree = sys.argv[17]
#If tree:
tree_type = str(sys.argv[18])
bootstrap_replications = str(sys.argv[19])
"""
#-----------------------------------------------------



def pro_link(query_proteins, hitlist_range, blast_database, smart_blast_, max_low_identity_seqs, min_low_identity_seqs, expected_min_identity, additional_hits, remove_gaps, cluster_seqs, similarity, smart_clustering_, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, align_seqs, generate_logo, weblogo_format, generate_tree, tree_type, bootstrap_replications):
    my_sequences = obtain_fasta_file(query_proteins)
    for seq_n, my_seq_record in enumerate(my_sequences):

        my_sequence_index = seq_n + 1

        os.mkdir("./outputs/protein_" + str(my_sequence_index))
        with open('./outputs/protein_' + str(my_sequence_index) + '/report.txt', 'w') as f:
            f.write('Protein_'+str(my_sequence_index)+" report:")
            f.writelines("Protein_id: " + str(my_seq_record.id))
        blast_filename = "./outputs/"+"protein_"+str(my_sequence_index)+"/blast_results.xml"
        found_sequences_fastafile= "./outputs/"+"protein_"+str(my_sequence_index) + "/found_sequences.fasta"

        if smart_blast_:
            s_blast(my_sequence_index, blast_database, hitlist_range, my_seq_record, blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits)
        else:
            blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
            parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity)

        cluster_results_fastafile = "./outputs/"+"protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + ".fasta"
        cluster_results_file = "./outputs/"+"protein_"+str(my_sequence_index) + "/cluster_results_"+ str(similarity)
        cluster_evaluation_file = "./outputs/"+"protein_"+str(my_sequence_index) +  "/cluster_results_evaluation_" + str(similarity)

        if cluster_seqs:
            if smart_clustering_:
                s_cluster(found_sequences_fastafile, my_seq_record, similarity, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            else:
                cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)


        if align_seqs:
            muscle_output= "./outputs/"+"protein_"+str(my_sequence_index) + "/cluster_results_evaluation_"  "aligned" ".fasta"
            align(cluster_results_fastafile, muscle_output)
            if generate_logo:
                weblogo_output = "./outputs/"+"protein_"+str(my_sequence_index) +'/logo'+'.'+str(weblogo_format)
                weblogo3(weblogo_format, muscle_output, weblogo_output)
            if generate_tree:
                mega_config_input = "./ProLink/modules/mega_configs/"+tree_type+"_"+bootstrap_replications+".mao"
                mega_output = "./outputs/"+"protein_"+str(my_sequence_index) +"/"+tree_type +"_"+bootstrap_replications+"_tree" 
                tree(mega_config_input, muscle_output, mega_output)

        else:
            print("Process finished (no clustering)")

pro_link(query_proteins, hitlist_range, blast_database, smart_blast_, max_low_identity_seqs, min_low_identity_seqs, expected_min_identity, additional_hits, remove_gaps, cluster_seqs, similarity, smart_clustering_, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, align_seqs, generate_logo, weblogo_format, generate_tree, tree_type, bootstrap_replications)
