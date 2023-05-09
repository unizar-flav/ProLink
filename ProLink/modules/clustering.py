
import logging
import os
import subprocess
from copy import copy
from textwrap import dedent

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def cluster(found_sequences_fastafile:str,
            similarity:float,
            cluster_results_file:str,
            cluster_evaluation_file:str,
            cluster_results_fastafile:str) -> int:
    '''
    Cluster sequences with ALFATClust

    Parameters
    ----------
    found_sequences_fastafile : str
        Path of the file with the sequences to cluster (FASTA format)
    similarity : float
        Similarity threshold for clustering
    cluster_results_file : str
        Path of the file with the clustering results
    cluster_evaluation_file : str
        Path of the file with the clustering evaluation
    cluster_results_fastafile : str
        Path of the file with the clustered sequences (FASTA format)

    Returns
    -------
    int
        Number of clusters
    '''
    logging.info(dedent(f"""
        -- Clustering sequences with ALFATClust
        Similarity:               {similarity}
        Input file:               '{found_sequences_fastafile}'
        Results file:             '{cluster_results_file}'
        Evaluation file:          '{cluster_evaluation_file}'
        """))
    alfatclust_cmd = ['alfatclust.py', '-i', found_sequences_fastafile, '-o', cluster_results_file, '-e', cluster_evaluation_file, '-ea' ,'-l' , str(similarity)]
    logging.debug(f"ALFATClust command: {alfatclust_cmd}")
    alfatclust_run = subprocess.run(alfatclust_cmd)
    if alfatclust_run.returncode != 0:
        logger.error(f"ERROR: ALFATClust failed")
        raise RuntimeError(f"ALFATClust failed")

    logger.info(f"Parsing clustering results")
    found_sequences = list(SeqIO.parse(found_sequences_fastafile, "fasta"))
    with open(cluster_results_file, 'r') as f:
        cluster_results = f.readlines()
    cluster_center_seqs = [cluster_results[n+1].strip() for n, line in enumerate(cluster_results) if line.startswith("#Cluster")]
    clustered_sequences = []
    for n_cluster, cluster_center_seq in enumerate(cluster_center_seqs, 1):
        logging.debug(f"Cluster {n_cluster}: {cluster_center_seq}")
        n_cluster_str = f"---C{n_cluster}"
        for seq_record in found_sequences:
            if seq_record.description.strip() == cluster_center_seq:
                seq_record_central = copy(seq_record)
                seq_record_central.id = cluster_center_seq + n_cluster_str
                seq_record_central.description = ""
                clustered_sequences.append(seq_record_central)
                break
        else:
            logger.warning(f"WARNING: Cluster {n_cluster} central sequence not found in the sequences file: {cluster_center_seq}")

    logger.debug(f"Writing clustered sequences to '{cluster_results_fastafile}'")
    SeqIO.write(clustered_sequences, cluster_results_fastafile, "fasta")
    n_clusters = len(clustered_sequences)
    logger.info(f"Clustering done. Number of clusters: {n_clusters}\n")
    return n_clusters

def p_cluster(found_sequences_fastafile:str,
              similarity:float,
              cluster_results_file:str,
              cluster_evaluation_file:str,
              cluster_results_fastafile:str,
              n_clusters_range:list[int]) -> int:
    '''
    Pro Cluster sequences with ALFATClust

    Parameters
    ----------
    found_sequences_fastafile : str
        Path of the file with the sequences to cluster (FASTA format)
    similarity : float
        Similarity threshold for clustering
    cluster_results_file : str
        Path of the file with the clustering results
    cluster_evaluation_file : str
        Path of the file with the clustering evaluation
    cluster_results_fastafile : str
        Path of the file with the clustered sequences (FASTA format)
    n_clusters_range : list[int]
        Range of number of clusters to cluster the sequences

    Returns
    -------
    int
        Number of clusters
    '''
    max_iter = 100
    for iteration in range(max_iter):
        logger.info(f"Pro Clustering iteration {iteration+1}\n")
        n_clusters = cluster(found_sequences_fastafile, similarity, cluster_results_file, cluster_evaluation_file, cluster_results_fastafile)
        if n_clusters_range[0] <= n_clusters <= n_clusters_range[1]:
            break
        logging.info(f"Number of clusters {n_clusters} not in range {n_clusters_range}")
        sign = 1 if n_clusters < n_clusters_range[0] else -1
        similarity += sign * 0.1
        logger.info(f"Clustering again with similarity = {similarity}")
        logger.debug(f"Removing files: '{cluster_results_file}', '{cluster_evaluation_file}', '{cluster_results_fastafile}'")
        os.remove(cluster_results_file)
        os.remove(cluster_evaluation_file)
        os.remove(cluster_results_fastafile)
    else:
        logger.error(f"ERROR: Pro Clustering failed: maximum number of iterations reached")
        raise Exception("Maximum number of iterations reached")
    logger.info(f"Pro Clustering done succesfully: {n_clusters} clusters found\n")
    return n_clusters
