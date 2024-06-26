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
from datetime import datetime, timezone

from . import __version__, ProLink_path, parameters_default
from .modules.blast import blast, blast_parse, p_blast
from .modules.clustering import cluster_mmseqs, cluster_pro
from .modules.obtaining_sequences import check_seq_in, get_seq
from .modules.pfam import fasta_to_dfasta
from .modules.subprocess_functions import align, tree
from .modules.trim import trim_align
from .modules.weblogo import weblogo3


logger = logging.getLogger()

def pro_link(query:str, parameters_default:dict = parameters_default, **parameters) -> None:
    '''
    Main function to run ProLink with a single query

    Parameters
    ----------
    query : str
        Sequence code of the protein to query
    parameters_default : dict, optional
        Default parameters for ProLink (def: taken from 'parameters_default')
    **parameters : dict
        Extra parameters to pass to ProLink
    '''

    # Add logger file handler if not present
    if not any(isinstance(handler, logging.FileHandler) for handler in logger.handlers):
        logger.addHandler(logging.FileHandler(f"{query}.log", mode='w'))

    time_start = datetime.now(timezone.utc)
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
    min_seq_id = float(parameters['min_seq_id'])
    pro_clustering_ = bool(parameters['pro_clustering_'])
    min_seq_id_step = float(parameters['min_seq_id_step'])
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
    output_dir = str(parameters['output_dir']) or f"{query}"

    # Manage output directory
    if os.path.exists(output_dir):
        logger.debug(f"Outputs directory already exists: {output_dir}. Overwriting.")
        os.removedirs(output_dir)
    logger.debug(f"Create outputs directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # Obtain sequence
    logger.info(f"Obtaining sequence for query protein: {query}")
    try:
        seq_record = get_seq(query, f"{output_dir}/my_sequence.fasta")[0]
        logger.info(f"\n> {seq_record.id} - {seq_record.description}")
        logger.info(f"{seq_record.seq}\n\n")
    except Exception as e:
        logger.debug("", exc_info=True)
        logger.error(f"ERROR: Obtaining sequences failed (Wrong query?): {e}")
        return

    # Process sequence
    try:
        blast_filename = f"{output_dir}/blast_results.xml"
        found_sequences_fastafile = f"{output_dir}/seqs_blast.fasta"
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

        check_seq_in(seq_record, found_sequences_fastafile, rewrite=True, spaces=False)

        if cluster_seqs:
            cluster_results = f"{output_dir}/seqs_cluster"
            cluster_results_fastafile = f"{cluster_results}.fasta"
            if pro_clustering_:
                logger.info(f"\n###  Pro Clustering  ###\n")
                cluster_pro(found_sequences_fastafile, cluster_results, [min_number_of_clusters_to_cluster_again, max_number_of_clusters_to_cluster_again], min_seq_id, min_seq_id_step)
            else:
                logger.info(f"\n###  Clustering  ###\n")
                cluster_mmseqs(found_sequences_fastafile, cluster_results, min_seq_id)
            sequences_fastafile = cluster_results_fastafile
            sequences_fastafile_pfam = f"{output_dir}/seqs_cluster_pfam.fasta"
            pfam_output = f"{output_dir}/seqs_cluster_pfam.txt"
            align_basename = f"{output_dir}/seqs_cluster_aligned"
        else:
            sequences_fastafile = found_sequences_fastafile
            sequences_fastafile_pfam = f"{output_dir}/seqs_blast_pfam.fasta"
            pfam_output = f"{output_dir}/seqs_blast_pfam.txt"
            align_basename = f"{output_dir}/seqs_blast_aligned"

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
                weblogo_output = f"{output_dir}/logo.{weblogo_format}"
                weblogo3(aligned_fastafile, weblogo_output, weblogo_format)
            if trim:
                logger.info("\nTrimming alignment")
                align_output_trim = f"{align_basename}_trim.fasta"
                trim_align(aligned_fastafile, align_output_trim)
                aligned_fastafile = align_output_trim
                if generate_logo:
                    logger.info("\nGenerating trimmed sequence logo")
                    weblogo_output_trim = f"{output_dir}/logo_trim.{weblogo_format}"
                    weblogo3(aligned_fastafile, weblogo_output_trim, weblogo_format)
            if generate_tree:
                logger.info("\nGenerating tree")
                mega_output = f"{aligned_fastafile}.mega"
                tree(tree_type, bootstrap_replications, aligned_fastafile, mega_output)
        else:
            logger.info("\nSkipping alignment (and logo and tree))")

    except Exception as e:
        logger.debug("", exc_info=True)
        logger.error(f"\nERROR: Fatal error on query {seq_record.id}. Aborting.")

    time_elapsed = datetime.now(timezone.utc) - time_start
    logger.info(f"\n\nProcess finished for query {query}\n" +
                f"Time elapsed: {time_elapsed.seconds//3600 + 24*time_elapsed.days}h {(time_elapsed.seconds//60)%60}m {time_elapsed.seconds%60}s\n\n")

    return

def pro_link_multiple(query_list:list[str], parameters_default:dict = parameters_default, **parameters) -> None:
    '''
    Run ProLink with multiple queries

    Parameters
    ----------
    query_list : list[str]
        List of sequence codes of the proteins to query
    parameters_default : dict, optional
        Default parameters for ProLink (def: taken from 'parameters_default')
    **parameters : dict
        Extra parameters to pass to ProLink
    '''

    # Add logger file handler if not present
    if not any(isinstance(handler, logging.FileHandler) for handler in logger.handlers):
        logger.addHandler(logging.FileHandler(f"ProLink.log", mode='w'))

    for query in query_list:
        try:
            pro_link(query, parameters_default, **parameters)
        except:
            continue
    return