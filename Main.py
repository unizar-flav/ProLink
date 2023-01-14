"""                                        
                                 __                                    
                                / /                                    
                               / /                                     
           ______ _____ _____ / /    __  ______ __ _                    
          / __  // ___// _  // /    / / / __  // ///                   
         / /_/ // /   / // // /___ / / / / / // _ \                     
        / ____//_/   /____//_____//_/ /_/ /_//_/ \_\                      
       / /                                                             
      / /                    Created by Víctor Sanz                                         
     /_/                    University of Zaragoza                                         
                                                                   
                                                                       
                                                                       
"""

from genericpath import isdir
import os
import subprocess
from .modules.obtaining_sequences import obtain_fasta_file
from .modules.blast import blast
from .modules.blast import parse
from .modules.blast import p_blast
from .modules.clustering import cluster
from .modules.clustering import p_cluster
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align
from .modules.subprocess_functions import weblogo3
from .modules.subprocess_functions import tree
import os
import sys

def pro_link(query_proteins, **parameters):
    
    hitlist_range = int(parameters['Blast']['hitlist_range'])
    blast_database = str(parameters['Blast']['blast_database'])
    pro_blast_ = bool(parameters['Blast']['pro_blast_'])
    max_low_identity_seqs = int(parameters['Blast']['max_low_identity_seqs'])
    min_low_identity_seqs = int(parameters['Blast']['min_low_identity_seqs'])
    expected_min_identity = float(parameters['Blast']['expected_min_identity'])
    additional_hits = int(parameters['Blast']['additional_hits'])
    remove_gaps = bool(parameters['Blast']['remove_gaps'])
    
    cluster_seqs = bool(parameters['Clustering']['cluster_seqs'])
    similarity = float(parameters['Clustering']['similarity'])
    pro_clustering_ = bool(parameters['Clustering']['pro_clustering_'])
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

        if pro_blast_:
            print("Pro BLAST")
            p_blast(my_sequence_index, blast_database, hitlist_range, my_seq_record, blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, outputs_dir)
        else:
            blast(hitlist_range, blast_database, blast_filename,  my_seq_record)
            parse(my_sequence_index, hitlist_range, blast_filename, remove_gaps, expected_min_identity, outputs_dir)

      
        if cluster_seqs:
            
            cluster_results_fastafile = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + ".fasta"
            cluster_results_file = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_"+ str(similarity)
            cluster_evaluation_file = "./" + outputs_dir +"/protein_"+str(my_sequence_index) +  "/cluster_results_evaluation_" + str(similarity)
            
            
            if pro_clustering_:
                print("Pro Clustering")
                p_cluster(found_sequences_fastafile, my_seq_record, similarity, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
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
                found_sequences_fastafile_pfam = "./" + outputs_dir + "/protein_"+str(my_sequence_index) + "/cluster_results_evaluation_" + str(similarity) + "_pfam.fasta"
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
                
                
          


