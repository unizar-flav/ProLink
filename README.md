# ***ProLink*** 

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1EhX_gO05Fmn_-ikvEkz90rf9S9p0imtp?usp=sharing)

## Overview
ProLink is a python script that allows to excecute multiple proteomic analysis tools automatically.

## Usage
The script is designed to be executed in Google Colab.

**Step 1:** Run the first cell in order to install the dependencies required for the script. It is only needed to do once everytime the execution environment is initialized. 

**Step 2:** Introduce the desired parameters in the form of the second cell and execute it.

***Parameters***

| Argument name                             | Description                                                                                      | 
| ----------------------------------------- | -------------------------------------------------------------------------------------------------|
| query_proteins                            | UniProt code of the query proteins. Eg: "ABQ62066.1, ABQ62091.1, ABQ62490.1".                    |
| hitlist_range                             | Number of found sequences obtained via Blast.                                                    |
| blast_database                            | Database used in blast.                                                                          |
| pro_blast_                                | Boolean parameter to select "Pro BLAST" or regular BLAST.                                        |
| cluster_seqs                              | Boolean parameter to select if clustering the sequences or not.                                  |
| similarity                                | Initial similarity treshold to group the sequences into clusters.                                |
| pro_clustering                            | Boolean parameter to select "pro clustering" or regular clustering.                              |
| check_pfam_domains                        | Boolean parameter to select if checking the Pfam domains of the sequences or not.                |
| align_seqs                                | Boolean parameter to select if aligning the sequences or not.                                    |
| generate_logo                             | Boolean parameter to select if generating a sequence logo or not.                                |
| generate_tree                             | Boolean parameter to select if generating a phylogenetic tree or not.                            |
| tree_type                                 | Type of phylogenetic tree. Needs to be either "NJ" (Neighbor joining) or "ML" (Maximum likehood).|


***Advanced parameters (in ProLink/parameters.cfg)***
| Argument name                             | Description                                                                                        | Default value|
| ----------------------------------------- | ---------------------------------------------------------------------------------------------------|--------------|
| max_low_identity_seqs                     | Maximum number of low identity seqs to find when using "Pro BLAST".                                |             1|
| min_low_identity_seqs                     | Minimum number of low identity seqs to find when using "Pro BLAST".                                |             1|
| expected_min_identity                     | Maximum identity percentage to consider a sequence a low identity seq when using "Pro BLAST".      |          0.25|
| additional_hits                           | Number of additional sequences to find when using "Pro BLAST".                                     |          2000|
| min_number_of_clusters_to_cluster_again   | Minimum number of clusters allowed when using "Pro clustering".                                    |           250|
| max_number_of_clusters_to_cluster_again   | Maximum number of clusters allowed when using "Pro clustering".                                    |           700|
| weblogo_format                            | Output format when using generate_logo.                                                            | png          |
| bootstrap_replications                    | Number of bootstrap replications when generating the tree. Needs to be 100, 250, 500, 1000 or 2000.|           100|
| outputs_dir                               | Name of the outputs directory.                                                                     | outputs      |

**Step 3:** Run the "Execute the script" cell. 

If the "download_outputs" parameter is marked, the outputs directory will be automatically downloaded to your device when finished.
As an additional option, one single file or directory can be downloaded using the "(Optional) Download outputs manually" cell.

## Advanced functions

### Pro BLAST

Pro BLAST allows to search for homologous sequences via the Blast tool, but making sure that low identity seqs are also represented in the output.

Firstly, a regular Blast with is launched, and hitlist_range seqs are obtained. 

If the number of low identity seqs in the found seqs are below the min_low_identity_seqs parameter, a new regular Blast will be launched but with hitlist_range += additional_hits.

Pro BLAST will stop when the number of low identity seqs reach max_low_identity_seqs, or when more than 10000 sequences are found (this limit is due to the computing capacity of the Google collab servers when clustering. In a future update, this value will be included as an advanced parameter).

### Pro clustering

ALFATClust is the tool that is used to cluster the sequences. It is a recent and relatively fast tool, but it does not let the user choose the number of clusters produced when grouping the sequences.

Pro clustering allows to obtain a number of clusters that is in a determined range. 

Firstly, a regular clustering with the determined similarity is executed. 

If the number of clusters is avobe the max_number_of_clusters_to_cluster_again value, the sequences will be clustered again but with similarity -= 0.1, in order to obtain an inferior number of clusters.

On the contrary, if the number of clusters is below the min_number_of_clusters_to_cluster_again, the sequences will be clustered again but with similarity += 0.1, in order to obtain a superior number of clusters.

Pro clustering will stop when the number of clusters is between min_number_of_clusters_to_cluster_again and max_number_of_clusters_to_cluster_again.

## Explicative diagram
![alt text](https://github.com/unizar-flav/ProLink/blob/main/diagram.png?raw=true)

## References

  > Chiu, J.K.H., Ong, R.TH. Clustering biological sequences with dynamic sequence similarity threshold. BMC Bioinformatics 23, 108 (2022).     https://doi.org/10.1186/s12859-022-04643-9
  
  > Ondov, B. D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology, 17(1), 132.
  
  > Steinegger, M. and J. Söding. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026.
  
  > Bakan A, Meireles LM, Bahar I. ProDy: Protein Dynamics Inferred from Theory and Experiments. Bioinformatics 2011 27(11):1575-1577.
  
  > Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 32(5), 1792-1797. https://doi.org/10.1093/nar/gkh340
  
  > Crooks GE, Hon G, Chandonia JM, Brenner SE WebLogo: A sequence logo generator, Genome Research, 14:1188-1190, (2004).
  
  > Koichiro Tamura, Glen Stecher, and Sudhir Kumar (2021) MEGA11: Molecular Evolutionary Genetics Analysis version 11. Molecular Biology and Evolution 38:3022-3027.
