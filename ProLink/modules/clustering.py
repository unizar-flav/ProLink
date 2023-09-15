
import logging
import os
import subprocess
from copy import copy
from tempfile import TemporaryDirectory
from textwrap import dedent

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def cluster_mmseqs(found_sequences_fastafile:str,
                   cluster_results:str,
                   min_seq_id:float=0.6) -> int:
    '''
    Cluster sequences with MMseqs2

    Parameters
    ----------
    found_sequences_fastafile : str
        Path of the file with the sequences to cluster (FASTA format)
    cluster_results : str
        Basename for the files with the clustering results (.tsv, .txt, .csv and .fasta)
    min_seq_id : float, optional
        Minimum sequence identity for clustering together (0-1) (def: 0.6)

    Returns
    -------
    int
        Number of clusters
    '''
    # min_seq_id between 0 and 1
    min_seq_id = min(max(min_seq_id, 0), 1)
    # output files
    cluster_results_tsv = f"{cluster_results}.tsv"
    cluster_results_txt = f"{cluster_results}.txt"
    cluster_results_csv = f"{cluster_results}.csv"
    cluster_results_fasta = f"{cluster_results}.fasta"
    # check if file exists and remove it
    for file in [cluster_results_tsv, cluster_results_txt, cluster_results_csv, cluster_results_fasta]:
        if os.path.exists(file):
            logger.debug(f"Removing file: '{file}'")
            os.remove(file)
    logging.info(dedent(f"""
        -- Clustering sequences with MMseqs2
        Minimum sequence identity:  {min_seq_id}
        Input file:                '{found_sequences_fastafile}'
        Results files:             '{cluster_results_txt}' / '{cluster_results_csv}' / '{cluster_results_fasta}'
        """))
    # perform the clustering
    with TemporaryDirectory() as tmp_dir:
        create_db_cmd = ['mmseqs', 'createdb', '--dbtype', '1', '--shuffle', '0', '--createdb-mode', '0', found_sequences_fastafile, os.path.join(tmp_dir, 'input_DB')]
        logging.debug(f"MMseqs2 create DB command: {' '.join(create_db_cmd)}")
        create_db_run = subprocess.run(create_db_cmd)
        cluster_cmd = ['mmseqs', 'cluster', '--cluster-mode', '0', '--min-seq-id', str(min_seq_id), os.path.join(tmp_dir, 'input_DB'), os.path.join(tmp_dir, 'cluster_DB'), os.path.join(tmp_dir, 'tmp')]
        logging.debug(f"MMseqs2 cluster command: {' '.join(cluster_cmd)}")
        cluster_run = subprocess.run(cluster_cmd)
        createtsv_cmd = ['mmseqs', 'createtsv', '--first-seq-as-repr', '0', '--full-header', '0', '--idx-seq-src', '0', os.path.join(tmp_dir, 'input_DB'), os.path.join(tmp_dir, 'input_DB'), os.path.join(tmp_dir, 'cluster_DB'), cluster_results_tsv]
        logging.debug(f"MMseqs2 createtsv command: {' '.join(createtsv_cmd)}")
        createtsv_run = subprocess.run(createtsv_cmd)
    # process the '.tsv' file
    found_sequences = list(SeqIO.parse(found_sequences_fastafile, "fasta"))
    with open(cluster_results_tsv, 'r') as f:
        cluster_results = f.readlines()
    clusters = {}
    with open(cluster_results_tsv, 'r') as f:
        for line in f:
            ref, seq = line.strip().split()
            clusters.setdefault(ref, [])
            for seq_record in found_sequences:
                if seq_record.id == seq:
                    clusters[ref].append(seq_record)
                    break
    logging.info(f"Clustering done. Number of clusters: {len(clusters)}")
    # write '.txt' file
    with open(cluster_results_txt, 'w') as f:
        for n_cluster, (ref, seqs) in enumerate(clusters.items(), 1):
            f.write(f"#Cluster {n_cluster}\n")
            for seq in seqs:
                f.write(f"{seq.description}\n")
    # write '.csv' file
    with open(cluster_results_csv, 'w') as f:
        f.write(f"# Cluster ID, No. of sequences, Center sequence\n")
        for n_cluster, (ref, seqs) in enumerate(clusters.items(), 1):
            f.write(f"{n_cluster}, {len(seqs)}, {seqs[0].description}\n")
    # write '.fasta' file
    clusters_fastafile = []
    for n_cluster, (ref, seqs) in enumerate(clusters.items(), 1):
        seq = copy(seqs[0])
        seq.id = f"{seq.description}---C{n_cluster}"
        seq.description = ""
        clusters_fastafile.append(seq)
    SeqIO.write(clusters_fastafile, cluster_results_fasta, "fasta")
    return len(clusters)

def cluster_alfatclust(found_sequences_fastafile:str,
                       similarity:float,
                       cluster_results_file:str,
                       cluster_evaluation_file:str,
                       cluster_results_fastafile:str) -> int:
    '''
    [Deprecated] Cluster sequences with ALFATClust

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
    raise DeprecationWarning("Clustering by ALFATClust is deprecated")
    logging.info(dedent(f"""
        -- Clustering sequences with ALFATClust
        Similarity:               {similarity}
        Input file:               '{found_sequences_fastafile}'
        Results file:             '{cluster_results_file}'
        Evaluation file:          '{cluster_evaluation_file}'
        """))
    alfatclust_cmd = ['alfatclust.py', '-i', found_sequences_fastafile, '-o', cluster_results_file, '-e', cluster_evaluation_file, '-ea' ,'-l' , str(similarity)]
    logging.debug(f"ALFATClust command: {' '.join(alfatclust_cmd)}")
    alfatclust_run = subprocess.run(alfatclust_cmd)
    if alfatclust_run.returncode != 0 or not os.path.exists(cluster_results_file):
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
              cluster_results:str,
              n_clusters_range:list[int],
              min_seq_id:float=0.6,
              min_seq_id_step:float=0.02) -> int:
    '''
    Pro Cluster sequences with MMseqs2

    Parameters
    ----------
    found_sequences_fastafile : str
        Path of the file with the sequences to cluster (FASTA format)
    cluster_results : str
        Basename for the files with the clustering results (.tsv, .txt, .csv and .fasta)
    n_clusters_range : list[int]
        Range of number of clusters to cluster the sequences
    min_seq_id : float, optional
        Minimum sequence identity for clustering together (0-1) (def: 0.6)
    min_seq_id_step : float, optional
        Step to increase or decrease the identity threshold (def: 0.02)

    Returns
    -------
    int
        Number of clusters
    '''
    #TODO: dinamic min_seq_id_step (guess from previous iteration results)
    found_sequences = list(SeqIO.parse(found_sequences_fastafile, "fasta"))
    n_clusters_range = sorted(n_clusters_range)
    n_clusters_range = [min(n_clusters_range[0], len(found_sequences)), min(n_clusters_range[1], len(found_sequences))]
    logger.info(f"Pro Clustering looking for {n_clusters_range[0]}-{n_clusters_range[1]} clusters")
    max_iter = 100
    min_seq_ids = set()
    for iteration in range(max_iter):
        logger.info(f"Pro Clustering iteration {iteration+1}\n")
        n_clusters = cluster_mmseqs(found_sequences_fastafile, cluster_results, min_seq_id)
        if n_clusters_range[0] <= n_clusters <= n_clusters_range[1]:
            break
        logging.info(f"Number of clusters {n_clusters} not in range {n_clusters_range}")
        if min_seq_id in min_seq_ids:
            logger.warning(f"WARNING: Pro Clustering failed converging (min_seq_id already used)")
            break
        sign = 1 if n_clusters < n_clusters_range[0] else -1
        min_seq_id += sign * min_seq_id_step
        min_seq_ids.add(min_seq_id)
        logger.info(f"Clustering again with min_seq_id = {min_seq_id}")
    else:
        logger.error(f"ERROR: Pro Clustering failed: maximum number of iterations reached")
        raise Exception("Maximum number of iterations reached")
    logger.info(f"Pro Clustering done succesfully: {n_clusters} clusters found\n")
    return n_clusters
