#!/usr/bin/env python3

import argparse

from .prolink import pro_link


def main():
    parser = argparse.ArgumentParser(
        prog='prolink',
        description='Excecute multiple proteomic analysis tools automatically')
    parser.add_argument('UNIPROT_ID', type=str, nargs='+',
        help='UniProt code of the protein(s) to query')
    args = parser.parse_args()
    
    query_proteins = ", ".join(args.UNIPROT_ID)
    pro_link(query_proteins)


if __name__ == '__main__':
    main()
