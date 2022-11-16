from genericpath import isdir
import os
import subprocess
from .modules.obtaining_sequences import obtain_fasta_file
from .modules.blast import blast
from .modules.blast import parse
from .modules.blast import s_blast
from .modules.clustering import cluster
from .modules.clustering import s_cluster
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align
from .modules.subprocess_functions import weblogo3
from .modules.subprocess_functions import tree
import os
import sys


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



def pro_link(query_proteins, **parameters):
    
    hitlist_range = int(parameters['Blast']['hitlist_range'])
    blast_database = str(parameters['Blast']['blast_database'])
    smart_blast_ = bool(parameters['Blast']['smart_blast_'])
    max_low_identity_seqs = int(parameters['Blast']['max_low_identity_seqs'])
    min_low_identity_seqs = int(parameters['Blast']['min_low_identity_seqs'])
    expected_min_identity = float(parameters['Blast']['expected_min_identity'])
    additional_hits = int(parameters['Blast']['additional_hits'])
    remove_gaps = bool(parameters['Blast']['remove_gaps'])
    
    cluster_seqs = bool(parameters['Clustering']['cluster_seqs'])
    similarity = float(parameters['Clustering']['similarity'])
    smart_clustering_ = bool(parameters['Clustering']['smart_clustering_'])
    min_number_of_clusters_to_cluster_again = int(parameters['Clustering']['min_number_of_clusters_to_cluster_again'])
    max_number_of_clusters_to_cluster_again = int(parameters['Clustering']['max_number_of_clusters_to_cluster_again'])
    
    check_pfam_domains = bool(parameters['Pfam domains']['check_pfam_domains'])
    
    align_seqs = bool(parameters['Alignment']['align_seqs'])
    
    generate_logo = bool(parameters['Weblogo']['generate_logo'])
    weblogo_format = str(parameters['Weblogo']['weblogo_format'])
    
    generate_tree = bool(parameters['Tree']['generate_tree'])
    tree_type = str(parameters['Tree']['tree_type'])
    bootstrap_replications = str(parameters['Tree']['bootstrap_replications'])
    
    outputs_dir = str(parameters['Outputs']['outputs_dir'])
    os.mkdir("./" + outputs_dir)

    my_sequences = obtain_fasta_file(query_proteins, outputs_dir)
    for seq_n, my_seq_record in enumerate(my_sequences):

        my_sequence_index = seq_n + 1

        os.mkdir("./"+ outputs_dir +"/protein_" + str(my_sequence_index))
        blast_filename = "./" + outputs_dir +"/protein_"+str(my_sequence_index)+"/blast_results.xml"
        found_sequences_fastafile= "./" + outputs_dir +"/protein_"+str(my_sequence_index) + "/found_sequences.fasta"

        if smart_blast_:
            print("Smart Blast")
            s_blast(my_sequence_index, blast_database, hitlist_range, my_seq_record, blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, outputs_dir)
        else:
            blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
            parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, outputs_dir)

      
        if cluster_seqs:
            
            cluster_results_fastafile = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + ".fasta"
            cluster_results_file = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_"+ str(similarity)
            cluster_evaluation_file = "./" + outputs_dir +"/protein_"+str(my_sequence_index) +  "/cluster_results_evaluation_" + str(similarity)
            
            
            if smart_clustering_:
                print("Smart Clustering")
                s_cluster(found_sequences_fastafile, my_seq_record, similarity, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            else:
                cluster(found_sequences_fastafile, my_seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            
            if check_pfam_domains:
                cluster_results_fastafile_pfam = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + "_pfam.fasta"
                fasta_to_dfasta(my_seq_record, cluster_results_fastafile, cluster_results_fastafile_pfam)
                cluster_results_fastafile = cluster_results_fastafile_pfam
            
            if align_seqs:
                muscle_output= "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_"  "aligned" ".fasta"
                print("Aligning sequences")
                align(cluster_results_fastafile, muscle_output)
                if generate_logo:
                    weblogo_output = "./" + outputs_dir + "/protein_"+str(my_sequence_index) +'/logo'+'.'+str(weblogo_format)
                    print("Generating sequence logo")
                    weblogo3(weblogo_format, muscle_output, weblogo_output)
                if generate_tree:
                    mega_config_input = "./ProLink/modules/mega_configs/"+tree_type+"_"+bootstrap_replications+".mao"
                    mega_output = "./" + outputs_dir + "/protein_"+str(my_sequence_index) +"/"+tree_type +"_"+bootstrap_replications+"_tree" 
                    print("Generating tree")
                    tree(mega_config_input, muscle_output, mega_output)

            else:
                print("Process finished (no alignment)")
        else:
            if check_pfam_domains:
                found_sequences_fastafile_pfam = = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + "_pfam.fasta"
                fasta_to_dfasta(my_seq_record, found_sequences_fastafile, found_sequences_fastafile_pfam)
                found_sequences_fastafile = found_sequences_fastafile_pfam
            
            if align_seqs:
                muscle_output= "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/found_sequences_"  "aligned" ".fasta"
                print("Aligning sequences")
                align(found_sequences_fastafile, muscle_output)
                if generate_logo:
                    weblogo_output = "./" + outputs_dir + "/protein_"+str(my_sequence_index) +'/logo'+'.'+str(weblogo_format)
                    print("Generating sequence logo")
                    weblogo3(weblogo_format, muscle_output, weblogo_output)
                if generate_tree:
                    mega_config_input = "./ProLink/modules/mega_configs/"+tree_type+"_"+bootstrap_replications+".mao"
                    mega_output = "./" + outputs_dir + "/protein_"+str(my_sequence_index) +"/"+tree_type +"_"+bootstrap_replications+"_tree" 
                    print("Generating tree")
                    tree(mega_config_input, muscle_output, mega_output)

            else:
                print("Process finished (no alignment)")    
                
    return outputs_dir
                
                
          


