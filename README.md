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
| query_proteins                            | UniProt code of the query proteins.                        |
                                              Eg: "ABQ62066.1, ABQ62091.1, ABQ62490.1"                   |
| hitlist_database                          | Desired minimum number of found sequences.                 |
