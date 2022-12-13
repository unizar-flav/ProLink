#!/usr/bin/env python3

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

import configparser
import os
from copy import deepcopy

from .modules.blast import blast, parse, s_blast
from .modules.clustering import cluster, s_cluster
from .modules.obtaining_sequences import obtain_fasta_file
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align, tree, weblogo3


ProLink_path = os.path.dirname(os.path.realpath(__file__))

# read default parameters from configuration file (as a plain dictionary)
config = configparser.ConfigParser()
config.read(os.path.join(ProLink_path, 'parameters.cfg'))
parameters_default = {}
for section in config.sections():
    parameters_default.update(dict(config[section]))


def pro_link(query_proteins: str, parameters_default: dict = parameters_default, **parameters):

    # assign default parameters that are not specified
    parameters_default = deepcopy(parameters_default)
    parameters_default.update(parameters)
    parameters = parameters_default

    # Blast
    hitlist_range = int(parameters['hitlist_range'])
    blast_database = str(parameters['blast_database'])
    smart_blast_ = bool(parameters['smart_blast_'])
    max_low_identity_seqs = int(parameters['max_low_identity_seqs'])
    min_low_identity_seqs = int(parameters['min_low_identity_seqs'])
    expected_min_identity = float(parameters['expected_min_identity'])
    additional_hits = int(parameters['additional_hits'])
    remove_gaps = bool(parameters['remove_gaps'])
    # Clustering
    cluster_seqs = bool(parameters['cluster_seqs'])
    similarity = float(parameters['similarity'])
    smart_clustering_ = bool(parameters['smart_clustering_'])
    min_number_of_clusters_to_cluster_again = int(parameters['min_number_of_clusters_to_cluster_again'])
    max_number_of_clusters_to_cluster_again = int(parameters['max_number_of_clusters_to_cluster_again'])
    # Pfam domains
    check_pfam_domains = bool(parameters['check_pfam_domains'])
    # Alignment
    align_seqs = bool(parameters['align_seqs'])
    # Weblogo
    generate_logo = bool(parameters['generate_logo'])
    weblogo_format = str(parameters['weblogo_format'])
    # Tree
    generate_tree = bool(parameters['generate_tree'])
    tree_type = str(parameters['tree_type'])
    bootstrap_replications = str(parameters['bootstrap_replications'])
    # Output
    outputs_dir = str(parameters['outputs_dir'])


    my_sequences = obtain_fasta_file(query_proteins, outputs_dir)
    for seq_n, my_seq_record in enumerate(my_sequences, 1):

        output_dir_n = f"./{outputs_dir}/protein_{seq_n}"
        os.makedirs(output_dir_n, exist_ok=True)

        blast_filename = f"{output_dir_n}/blast_results.xml"
        found_sequences_fastafile = f"{output_dir_n}/found_sequences.fasta"

        if smart_blast_:
            print("Smart Blast")
            s_blast(seq_n, blast_database, hitlist_range, my_seq_record, blast_filename, found_sequences_fastafile,
                    remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, outputs_dir)
        else:
            blast(hitlist_range, blast_database,
                  blast_filename,  my_seq_record)
            parse(seq_n, hitlist_range, blast_filename,
                  remove_gaps, expected_min_identity, outputs_dir)

        if cluster_seqs:

            cluster_results_file = f"{output_dir_n}/cluster_results_{similarity}"
            cluster_results_fastafile = f"{output_dir_n}/cluster_results_evaluation_{similarity}.fasta"
            cluster_evaluation_file = f"{output_dir_n}/cluster_results_evaluation_{similarity}"

            if smart_clustering_:
                print("Smart Clustering")
                s_cluster(found_sequences_fastafile, my_seq_record, similarity, min_number_of_clusters_to_cluster_again,
                          max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            else:
                cluster(found_sequences_fastafile, my_seq_record, similarity,
                        cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)

            if check_pfam_domains:
                cluster_results_fastafile_pfam = f"{output_dir_n}/cluster_results_evaluation_{similarity}_pfam.fasta"
                fasta_to_dfasta(my_seq_record, cluster_results_fastafile, cluster_results_fastafile_pfam)
                cluster_results_fastafile = cluster_results_fastafile_pfam

            if align_seqs:
                print("Aligning sequences")
                muscle_output = f"{output_dir_n}/cluster_results_evaluation_{similarity}_aligned.fasta"
                align(cluster_results_fastafile, muscle_output)
                if generate_logo:
                    print("Generating sequence logo")
                    weblogo_output = f"{output_dir_n}/logo.{weblogo_format}"
                    weblogo3(weblogo_format, muscle_output, weblogo_output)
                if generate_tree:
                    print("Generating tree")
                    mega_config_input = f"{ProLink_path}/modules/mega_configs/{tree_type}_{bootstrap_replications}.mao"
                    mega_output = f"{output_dir_n}/{tree_type}_{bootstrap_replications}_tree"
                    tree(mega_config_input, muscle_output, mega_output)
            else:
                print("Process finished (no alignment)")

        else:
            if check_pfam_domains:
                found_sequences_fastafile_pfam = f"{output_dir_n}/cluster_results_evaluation_{similarity}_pfam.fasta"
                fasta_to_dfasta(my_seq_record, found_sequences_fastafile, found_sequences_fastafile_pfam)
                found_sequences_fastafile = found_sequences_fastafile_pfam

            if align_seqs:
                print("Aligning sequences")
                muscle_output = f"{output_dir_n}/found_sequences_aligned.fasta"
                align(found_sequences_fastafile, muscle_output)
                if generate_logo:
                    print("Generating sequence logo")
                    weblogo_output = f"{output_dir_n}/logo.{weblogo_format}"
                    weblogo3(weblogo_format, muscle_output, weblogo_output)
                if generate_tree:
                    print("Generating tree")
                    mega_config_input = f"{ProLink_path}/modules/mega_configs/{tree_type}_{bootstrap_replications}.mao"
                    mega_output = f"{output_dir_n}/{tree_type}_{bootstrap_replications}_tree"
                    tree(mega_config_input, muscle_output, mega_output)
            else:
                print("Process finished (no alignment)")

    return outputs_dir
