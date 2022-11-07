# ProLink
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1EhX_gO05Fmn_-ikvEkz90rf9S9p0imtp?usp=sharing)

## Overview
ProLink is a python script that allows to excecute multiple proteomic analysis tools automatically.

## Usage
The script is designed to be executed in Google Colab.

**Step 1:** Run the first cell in order to install the environment required for the script. It is only neededed to do once everytime the execution environment is initialized. 

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
| min_low_identity_seqs                     | Minimum number of low identity seqs to find when using "smart blast".|
