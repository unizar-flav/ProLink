#!/usr/bin/env python3

r"""
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

import logging
import os
from copy import deepcopy
from datetime import datetime

from . import __version__, ProLink_path, parameters_default
from .modules.blast import blast, blast_parse, p_blast
from .modules.clustering import cluster, p_cluster
from .modules.obtaining_sequences import check_seq_in, get_seq
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align, tree
from .modules.trim import trim_align
from .modules.weblogo import weblogo3


logger = logging.getLogger()

def pro_link(query_proteins:list[str], parameters_default:dict = parameters_default, **parameters) -> None:

    time_start = datetime.utcnow()
    logger.info(f"ProLink v{__version__} started at {time_start.strftime('%Y-%m-%d %H:%M:%S')} UTC\n")

    logger.debug(f"ProLink path: {ProLink_path}")
    logger.debug(f"Default parameters: {parameters_default}")

    # assign default parameters that are not specified
    parameters_default = deepcopy(parameters_default)
    parameters_default.update(parameters)
    parameters = parameters_default
    logger.debug(f"Parameters: {parameters}")

    # Blast
    hitlist_size = int(parameters['hitlist_size'])
    blast_database = str(parameters['blast_database'])
    blast_local = bool(parameters['blast_local'])
    pro_blast_ = bool(parameters['pro_blast_'])
    length_restrict = bool(parameters['length_restrict'])
    length_margin = float(parameters['length_margin'])
    max_low_identity_seqs = int(parameters['max_low_identity_seqs'])
    min_low_identity_seqs = int(parameters['min_low_identity_seqs'])
    expected_min_identity = float(parameters['expected_min_identity'])
    additional_hits = int(parameters['additional_hits'])
    # Clustering
    cluster_seqs = bool(parameters['cluster_seqs'])
    similarity = float(parameters['similarity'])
    pro_clustering_ = bool(parameters['pro_clustering_'])
    similarity_step = float(parameters['similarity_step'])
    min_number_of_clusters_to_cluster_again = int(parameters['min_number_of_clusters_to_cluster_again'])
    max_number_of_clusters_to_cluster_again = int(parameters['max_number_of_clusters_to_cluster_again'])
    # Pfam domains
    check_pfam_domains = bool(parameters['check_pfam_domains'])
    # Alignment
    align_seqs = bool(parameters['align_seqs'])
    trim = bool(parameters['trim'])
    # Weblogo
    generate_logo = bool(parameters['generate_logo'])
    weblogo_format = str(parameters['weblogo_format'])
    # Tree
    generate_tree = bool(parameters['generate_tree'])
    tree_type = str(parameters['tree_type'])
    bootstrap_replications = int(parameters['bootstrap_replications'])
    # Output
    outputs_dir = str(parameters['outputs_dir'])

    logger.debug(f"Create outputs directory: {outputs_dir}")
    os.makedirs(outputs_dir, exist_ok=True)
    logger.info(f"Obtaining sequences for query proteins:  {', '.join(query_proteins)}")
    try:
        my_sequences = get_seq(query_proteins, f"{outputs_dir}/my_sequences.fasta")
        for seq_record in my_sequences:
            logger.info(f"\n> {seq_record.id} - {seq_record.description}")
            logger.debug(f"{seq_record.seq}")
    except Exception as e:
        logger.debug("", exc_info=True)
        logger.error(f"ERROR: Obtaining sequences failed (Wrong query?): {e}")
        return

    for seq_n, seq_record in enumerate(my_sequences, 1):

        try:
            logger.info(f"\n\n@@@@  Sequence {seq_n:02d}: {seq_record.id}\n")

            output_dir_n = f"./{outputs_dir}/{seq_n:02d}_{seq_record.id}"
            logger.debug(f"Create output directory for sequence {seq_n:02d}: {output_dir_n}")
            os.makedirs(output_dir_n, exist_ok=True)

            blast_filename = f"{output_dir_n}/blast_results.xml"
            found_sequences_fastafile = f"{output_dir_n}/found_sequences.fasta"
            if length_restrict:
                length_margin_seq = int(length_margin*len(seq_record.seq))
                length_range = [len(seq_record.seq) - length_margin_seq, len(seq_record.seq) + length_margin_seq]
                logger.debug(f"Length restriction range: {length_range[0]} - {length_range[1]}")
            else:
                length_range = []

            if pro_blast_:
                logger.info(f"\n###  Pro BLAST  ###\n")
                p_blast(seq_record, blast_filename, found_sequences_fastafile, expected_min_identity, min_low_identity_seqs, max_low_identity_seqs, additional_hits, hitlist_size, length_range, blast_database, blast_local)
            else:
                logger.info(f"\n###  BLAST  ###\n")
                blast(seq_record, blast_filename, blast_database, hitlist_size, blast_local)
                blast_parse(blast_filename, found_sequences_fastafile, expected_min_identity, True, max_low_identity_seqs, None, length_range)

            check_seq_in(seq_record, found_sequences_fastafile)

            if cluster_seqs:
                cluster_results_file = f"{output_dir_n}/cluster_seqs.txt"
                cluster_evaluation_file = f"{output_dir_n}/cluster_seqs_eval.csv"
                cluster_results_fastafile = f"{output_dir_n}/cluster_seqs_eval.fasta"
                if pro_clustering_:
                    logger.info(f"\n###  Pro Clustering  ###\n")
                    p_cluster(found_sequences_fastafile, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile, [min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again], similarity_step)
                else:
                    logger.info(f"\n###  Clustering  ###\n")
                    cluster(found_sequences_fastafile, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
                sequences_fastafile = cluster_results_fastafile
                sequences_fastafile_pfam = f"{output_dir_n}/cluster_seqs_eval_pfam.fasta"
                pfam_output = f"{output_dir_n}/cluster_seqs_eval_pfam.txt"
                align_basename = f"{output_dir_n}/cluster_seqs_eval_aligned"
            else:
                sequences_fastafile = found_sequences_fastafile
                sequences_fastafile_pfam = f"{output_dir_n}/found_seqs_pfam.fasta"
                pfam_output = f"{output_dir_n}/found_seqs_pfam.txt"
                align_basename = f"{output_dir_n}/found_seqs_aligned"

            if check_pfam_domains:
                logger.info("\nChecking Pfam domains")
                try:
                    fasta_to_dfasta(seq_record, sequences_fastafile, sequences_fastafile_pfam, pfam_output)
                    sequences_fastafile = sequences_fastafile_pfam
                except:
                    logger.debug("", exc_info=True)
                    logger.warning("WARNING: Errors while checking Pfam domains. Skipping.")

            if align_seqs:
                logger.info("\nAligning sequences")
                aligned_fastafile = f"{align_basename}.fasta"
                align(sequences_fastafile, aligned_fastafile)
                if generate_logo:
                    logger.info("\nGenerating sequence logo")
                    weblogo_output = f"{output_dir_n}/logo.{weblogo_format}"
                    weblogo3(aligned_fastafile, weblogo_output, weblogo_format)
                if trim:
                    logger.info("\nTrimming alignment")
                    align_output_trim = f"{align_basename}_trim.fasta"
                    trim_align(aligned_fastafile, align_output_trim)
                    aligned_fastafile = align_output_trim
                    if generate_logo:
                        logger.info("\nGenerating trimmed sequence logo")
                        weblogo_output_trim = f"{output_dir_n}/logo_trim.{weblogo_format}"
                        weblogo3(aligned_fastafile, weblogo_output_trim, weblogo_format)
                if generate_tree:
                    logger.info("\nGenerating tree")
                    mega_output = f"{aligned_fastafile}.mega"
                    tree(tree_type, bootstrap_replications, aligned_fastafile, mega_output)
            else:
                logger.info("\nSkipping alignment (and logo and tree))")

        except Exception as e:
            logger.debug("", exc_info=True)
            logger.error(f"ERROR: Fatal error on query {seq_record.id}. Aborting.")
            continue

        logger.info(f"\nProcess finished for query {seq_record.id}\n\n")

    time_elapsed = datetime.utcnow() - time_start
    logger.info(f"End of ProLink. Process finished for all queries.\n" +
                f"Time elapsed: {time_elapsed.seconds//3600 + 24*time_elapsed.days}h {(time_elapsed.seconds//60)%60}m {time_elapsed.seconds%60}s\n\n")
