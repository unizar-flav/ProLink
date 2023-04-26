#!/usr/bin/env python3

import argparse

from .prolink import pro_link


def main():
    parser = argparse.ArgumentParser(
        prog='prolink',
        description='Excecute multiple proteomic analysis tools automatically')
    parser.add_argument('UNIPROT_ID', type=str, nargs='+',
        help='UniProt code of the protein(s) to query')
    parser.add_argument('--options', metavar="opt:val", nargs='+', type=str,
        help='Extra options to pass: option name and value, separated by a colon (e.g. --options option1:value1 option2:value2)')
    args = parser.parse_args()
    
    query_proteins = ", ".join(args.UNIPROT_ID)
    options = dict()
    if args.options:
        for opt in args.options:
            opt_val= opt.split(":")
            if len(opt_val) != 2:
                print(f"ERROR: Options must be specified as 'option:value'. Got '{opt}'")
                exit(1)
            options[opt_val[0]] = opt_val[1]
    pro_link(query_proteins, **options)


if __name__ == '__main__':
    main()
