# ***ProLink*** 

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1EhX_gO05Fmn_-ikvEkz90rf9S9p0imtp?usp=sharing)

## Overview
ProLink is a python script that allows to excecute multiple proteomic analysis tools automatically.

## Usage
The script is designed to be executed in Google Colab.

**Step 1:** Run the first cell in order to install the environment required for the script. It is only needed to do once everytime the execution environment is initialized. 

**Step 2:** Introduce the desired parameters in the form of the second cell and execute it.

***Parameters***

| Argument name                             | Description                                                |
| ----------------------------------------- | ---------------------------------------------------------- |
| query_proteins                            | UniProt code of the query proteins. Eg: "ABQ62066.1, ABQ62091.1, ABQ62490.1"|
| blast_database                            | Database used in blast. Eg: "refseq_protein".              |
| smart_blast_                              | Boolean parameter to select "smart blast" or regular blast.|
| max_low_identity_seqs                     | Maximum number of low identity seqs to find when using "smart blast".|
| min_low_identity_seqs                     | Minimum number of low identity seqs to find when using "smart blast".|
| expected_min_identity                     | Maximum identity percentage to consider a sequence a low identity seq when using "smart blast".|
| additional_hits                           | Number of additional sequences to find when using "smart blast".|
| cluster_seqs                              | Boolean parameter to select if clustering the sequences or not.|
| similarity                                | Initial similarity treshold to group the sequences into clusters.|
| smart_clustering                          | Boolean parameter to select "smart clustering" or regular clustering.|
| min_number_of_clusters_to_cluster_again   | Minimum number of clusters allowed when using "smart clustering".|
| max_number_of_clusters_to_cluster_again   | Maximum number of clusters allowed when using "smart clustering".|
| cluster_seqs                              | Boolean parameter to select if aligning the sequences or not.|
| generate_logo                             | Boolean parameter to select if generating a sequence logo or not.|
| weblogo_format                            | Output format when using generate_logo. Eg: "png".|
| generate_tree                             | Boolean parameter to select if generating a phylogenetic tree or not.|
| tree_type                                 | Type of phylogenetic tree. Needs to be either "NJ" (Neighbor joining) or "ML" (Maximum likehood).|
| weblogo_format                            | Number of bootstrap replications when generating the tree.|
