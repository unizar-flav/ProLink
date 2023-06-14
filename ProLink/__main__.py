#!/usr/bin/env python3

import argparse
import logging

import yaml

from .prolink import __version__, parameters_default, pro_link


logger = logging.getLogger()

def main():
    parser = argparse.ArgumentParser(
        prog='prolink',
        description='Excecute multiple proteomic analysis tools automatically')
    parser.add_argument('-v', '--version', action='version', version=f'ProLink v{__version__}')
    parser.add_argument('QUERY_CODE', type=str, nargs='+',
        help='Sequence codes of the protein(s) to query')
    parser.add_argument('-f', '--file', metavar='.yml', type=str,
        help='File in YAML format with extra options to pass')
    parser.add_argument('--opt', metavar="opt=val", nargs='+', type=str,
        help='Extra options to pass: option name and value, separated by an equal (e.g. --opt opt1=val1 opt2=val2)')
    parser.add_argument('-o', '--outputs_dir', metavar='<>', type=str,
        help=f'Directory to write the outputs (def: \'{parameters_default["outputs_dir"]}\'))')
    parser.add_argument('--verbose', action='store_true',
        help='Verbose mode')
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    options = dict()

    # options from file
    if args.file:
        with open(args.file, 'r') as f:
            options_yaml = yaml.safe_load(f)
        # flatten options if necessary
        for k, v in options_yaml.items():
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    options[k2] = v2
            else:
                options[k] = v

    # options from command line
    if args.opt:
        for opt in args.opt:
            opt_val = opt.split("=")
            if len(opt_val) != 2:
                logger.error(f"ERROR: Options must be specified as 'option=value'. Got '{opt}'")
                exit(1)
            if opt_val[0] in parameters_default and isinstance(parameters_default[opt_val[0]], bool):
                if opt_val[1].lower() in ('y', 'yes', 't', 'true', 'on'):
                    opt_val[1] = True
                elif opt_val[1].lower() in ('n', 'no', 'f', 'false', 'off'):
                    opt_val[1] = False
            options[opt_val[0]] = opt_val[1]

    # outputs directory
    if args.outputs_dir:
        options['outputs_dir'] = args.outputs_dir

    pro_link(args.QUERY_CODE, **options)


if __name__ == '__main__':
    main()
