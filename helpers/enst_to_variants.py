#!/usr/bin/env python3

import sys
import csv
import json
from argparse import ArgumentParser
from collections import defaultdict


def parse_cmd_arguments():
    """
    Parse command-line arguments.

    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument(
        '--enst', '-e', dest='enst', type=str, required=True,
        help='''Ensembl stable transcript ID.'''
    )
    parser.add_argument(
        '--database', '-d', dest='database', type=str, required=True,
        help='''Database of variants in JSON format.'''
    )
    parser.add_argument(
        '--prefix', '-p', dest='prefix', type=str, required=True,
        help='Prefix to be prepended to the output filenames.'
    )
    return parser.parse_args()


def main():
    """
    This is where the variant extraction happens.

    Returns
    -------

    """
    # parse command-line arguments
    cmd_args = parse_cmd_arguments()

    # parse variant database
    with open(cmd_args.database, 'rt') as ipf:
        variant_db = json.load(ipf)

    enst_id = cmd_args.enst.strip()
    try:
        enst_variants = variant_db[enst_id]['variants']
    except KeyError:
        print('Transcript {} does not exist in the variant database.'.format(enst_id))
        sys.exit(1)

    # count variants
    syn_variants = defaultdict(int)
    mis_variants = defaultdict(int)
    for variant_info in enst_variants:
        variant, _, _ = variant_info
        w = variant[0]
        v = variant[-1]
        pos = int(variant[1:-1])
        if w == v:
            syn_variants[pos] += 1
        else:
            mis_variants[pos] += 1

    # write variants to file
    with open(cmd_args.prefix + '_' + enst_id + '_syn.csv', 'wt') as opf:
        syn_writer = csv.writer(opf)
        syn_writer.writerows([[pos, c] for pos, c in syn_variants.items()])
    with open(cmd_args.prefix + '_' + enst_id + '_mis.csv', 'wt') as opf:
        mis_writer = csv.writer(opf)
        mis_writer.writerows([[pos, c] for pos, c in mis_variants.items()])


if __name__ == '__main__':
    main()
