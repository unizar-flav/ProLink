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

import os
from copy import deepcopy

from . import ProLink_path, parameters_default
from .modules.blast import blast, p_blast, parse
from .modules.clustering import cluster, p_cluster
from .modules.obtaining_sequences import get_seq
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align, tree, weblogo3


def pro_link(query_proteins:list[str], parameters_default:dict = parameters_default, **parameters):

    # assign default parameters that are not specified
    parameters_default = deepcopy(parameters_default)
    parameters_default.update(parameters)
    parameters = parameters_default

    # Blast
    hitlist_size = int(parameters['hitlist_size'])
    blast_database = str(parameters['blast_database'])
    pro_blast_ = bool(parameters['pro_blast_'])
    max_low_identity_seqs = int(parameters['max_low_identity_seqs'])
    min_low_identity_seqs = int(parameters['min_low_identity_seqs'])
    expected_min_identity = float(parameters['expected_min_identity'])
    additional_hits = int(parameters['additional_hits'])
    remove_gaps = bool(parameters['remove_gaps'])
    # Clustering
    cluster_seqs = bool(parameters['cluster_seqs'])
    similarity = float(parameters['similarity'])
    pro_clustering_ = bool(parameters['pro_clustering_'])
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

    os.makedirs(outputs_dir, exist_ok=True)
    my_sequences = get_seq(query_proteins, f"{outputs_dir}/my_sequences.fasta")

    for seq_n, seq_record in enumerate(my_sequences, 1):

        output_dir_n = f"./{outputs_dir}/protein_{seq_n}"
        os.makedirs(output_dir_n, exist_ok=True)

        blast_filename = f"{output_dir_n}/blast_results.xml"
        found_sequences_fastafile = f"{output_dir_n}/found_sequences.fasta"

        if pro_blast_:
            print("Pro BLAST")
            p_blast(seq_record, blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, hitlist_size, database=blast_database)
        else:
            blast(seq_record, blast_filename, hitlist_size=hitlist_size, database=blast_database)
            parse(blast_filename, found_sequences_fastafile, remove_gaps, expected_min_identity)

        if cluster_seqs:
            cluster_results_file = f"{output_dir_n}/cluster_results_{similarity}"
            cluster_results_fastafile = f"{output_dir_n}/cluster_results_evaluation_{similarity}.fasta"
            cluster_evaluation_file = f"{output_dir_n}/cluster_results_evaluation_{similarity}"
            if pro_clustering_:
                print("Pro Clustering")
                p_cluster(found_sequences_fastafile, seq_record, similarity, min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            else:
                cluster(found_sequences_fastafile, seq_record, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
            sequences_fastafile = cluster_results_fastafile
            sequences_fastafile_pfam = f"{output_dir_n}/cluster_results_evaluation_{similarity}_pfam.fasta"
            muscle_output = f"{output_dir_n}/cluster_results_evaluation_{similarity}_aligned.fasta"
        else:
            sequences_fastafile = found_sequences_fastafile
            sequences_fastafile_pfam = f"{output_dir_n}/found_sequences_pfam.fasta"
            muscle_output = f"{output_dir_n}/found_sequences_aligned.fasta"

        if check_pfam_domains:
            fasta_to_dfasta(seq_record, sequences_fastafile, sequences_fastafile_pfam)
            sequences_fastafile = sequences_fastafile_pfam

        if align_seqs:
            print("Aligning sequences")
            align(sequences_fastafile, muscle_output)
            if generate_logo:
                print("Generating sequence logo")
                weblogo_output = f"{output_dir_n}/logo.{weblogo_format}"
                weblogo3(weblogo_format, muscle_output, weblogo_output)
            if generate_tree:
                print("Generating tree")
                mega_output = f"{output_dir_n}/cluster_results_evaluation_{similarity}_aligned.mega"
                tree(tree_type, bootstrap_replications, muscle_output, mega_output)
        else:
            print("Skipping alignment")

        print("Process finished")
