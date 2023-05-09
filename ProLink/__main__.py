#!/usr/bin/env python3

import argparse
import logging

from .prolink import pro_link


logger = logging.getLogger()

def main():
    parser = argparse.ArgumentParser(
        prog='prolink',
        description='Excecute multiple proteomic analysis tools automatically')
    parser.add_argument('UNIPROT_ID', type=str, nargs='+',
        help='UniProt code of the protein(s) to query')
    parser.add_argument('--options', metavar="opt=val", nargs='+', type=str,
        help='Extra options to pass: option name and value, separated by an equal (e.g. --options option1=value1 option2=value2)')
    parser.add_argument('--verbose', action='store_true',
        help='Verbose mode')
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    options = dict()
    if args.options:
        for opt in args.options:
            opt_val = opt.split("=")
            if len(opt_val) != 2:
                logger.error(f"ERROR: Options must be specified as 'option=value'. Got '{opt}'")
                exit(1)
            if opt_val[1].lower() in {"f", "false", "no", "off", "0"}:
                opt_val[1] = False
            options[opt_val[0]] = opt_val[1]

    pro_link(args.UNIPROT_ID, **options)


if __name__ == '__main__':
    main()
