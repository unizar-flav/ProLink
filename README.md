# ***ProLink***

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/unizar-flav/ProLink/blob/master/ProLink.ipynb)


## Overview
ProLink is a python script that allows to execute multiple proteomic analysis tools automatically.


## Usage
The script is designed to be executed in Google Colab. Although it can be installed and executed locally.

**Step 1:** Run the first cell in order to install the dependencies required for the script. It is only needed to do once everytime the execution environment is initialized.

**Step 2:** Introduce the desired parameters in the form of the second cell and execute it.

***Parameters***
| Argument name                             | Description                                                                                      |
| ----------------------------------------- | -------------------------------------------------------------------------------------------------|
| query_proteins                            | Protein sequence code to query. Will be searched on [NCBI](www.ncbi.nlm.nih.gov/Entrez).         |
| hitlist_size                              | Number of found sequences to obtain after BLAST.                                                 |
| blast_database                            | Database used in BLAST.                                                                          |
| pro_blast_                                | Boolean to select [*Pro BLAST*](#pro-blast) or regular *BLAST*.                                  |
| length_restrict                           | Boolean to resctrict the length of the found sequences with respect the query.                   |
| length_margin                             | Number to multiply the query length and restrict the min and max length of the found sequences.  |
| cluster_seqs                              | Boolean to select if clustering the sequences.                                                   |
| similarity                                | Initial similarity treshold to group the sequences into clusters.                                |
| pro_clustering_                           | Boolean to select [*Pro Clustering*](#pro-clustering) or regular clustering.                     |
| similarity_step                           | Step to increase or decrease the similarity threshold while *Pro Clustering*.                    |
| min_number_of_clusters_to_cluster_again   | Minimum number of clusters allowed when using *Pro Clustering*.                                  |
| max_number_of_clusters_to_cluster_again   | Maximum number of clusters allowed when using *Pro Clustering*.                                  |
| check_pfam_domains                        | Boolean to select if checking the Pfam domains of the sequences.                                 |
| align_seqs                                | Boolean to select if aligning the sequences.                                                     |
| trim                                      | Boolean to select if trimming the alignment retaining phylogenetically-informative sites.        |
| generate_logo                             | Boolean to select if generating a sequence logo.                                                 |
| generate_tree                             | Boolean to select if generating a phylogenetic tree.                                             |
| tree_type                                 | Type of phylogenetic tree. Either "NJ" (Neighbor Joining) or "ML" (Maximum Likehood).            |

***Advanced parameters (in ProLink/parameters.yaml)***
| Argument name                             | Description                                                                                        | Default value |
| ----------------------------------------- | ---------------------------------------------------------------------------------------------------|---------------|
| max_low_identity_seqs                     | Maximum number of low identity seqs to find when using "Pro BLAST".                                |             1 |
| min_low_identity_seqs                     | Minimum number of low identity seqs to find when using "Pro BLAST".                                |             1 |
| expected_min_identity                     | Maximum identity percentage to consider a sequence a low identity seq when using "Pro BLAST".      |          0.25 |
| additional_hits                           | Number of additional sequences to find when using "Pro BLAST".                                     |          2000 |
| weblogo_format                            | Output format when using generate_logo.                                                            |         'png' |
| bootstrap_replications                    | Number of bootstrap replications when generating the tree. Needs to be 100, 250, 500, 1000 or 2000.|           100 |
| outputs_dir                               | Name of the outputs directory.                                                                     |     'outputs' |

**Step 3:** Run the "Execute the script" cell. This may take a while to be completed. Check the box `extra_verbose` for an enriched output beyond the standard.

**Step 4:** Download the results as a zip file.


## Advanced functions

### Pro BLAST

*Pro BLAST* allows to search for homologous sequences via the BLAST tool, but making sure that low identity seqs are also represented in the output.

Firstly, a regular BLAST is launched and `hitlist_size` seqs are obtained.

If the number of low identity seqs in the found seqs are below the `min_low_identity_seqs` parameter, a new regular BLAST will be launched but with hitlist_size += `additional_hits`.

It will stop when the number of low identity seqs reach `max_low_identity_seqs`, or when more than 10000 sequences are found (this limit is due to the computing capacity of the Google collab servers when clustering. In a future update, this value will be included as an advanced parameter).

### Pro clustering

ALFATClust is the tool that is used to cluster the sequences. It is a recent and relatively fast tool, but it does not let the user choose the number of clusters produced when grouping the sequences.

*Pro Clustering* allows to obtain a number of clusters in a determined range.

Firstly, a regular clustering with the determined similarity is executed.

If the number of clusters is avobe the `max_number_of_clusters_to_cluster_again` value, the sequences will be clustered again but with similarity -= `similarity_step`, in order to obtain an inferior number of clusters. On the contrary, if the number of clusters is below the `min_number_of_clusters_to_cluster_again`, the similarity requested will be increased.


## References

  > Chiu, J.K.H., Ong, R.TH. (2022). Clustering biological sequences with dynamic sequence similarity threshold. BMC Bioinformatics 23, 108.

  > Ondov, B. D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132.

  > Steinegger, M., Söding., J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026.

  > Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 32(5), 1792-1797.

  > Crooks, G. E., et al. (2004).  WebLogo: A sequence logo generator, Genome Research, 14:1188-1190.

  > Tamura, K., et al. (2021) MEGA11: Molecular Evolutionary Genetics Analysis version 11. Molecular Biology and Evolution 38:3022-3027.
