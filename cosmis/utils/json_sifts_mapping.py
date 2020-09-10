#!/usr/bin/env python3


import os
import json
import gzip
from collections import defaultdict
from argparse import ArgumentParser


def main():
    """
    """
    # parse command-line arguments
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
            help='Input file containing the mappings.')
    parser.add_argument('-o', '--output', dest='output', required=False, type=str,
            help='Output JSON file containing the mappings.')
    args = parser.parse_args()

    #
    uniprot_pdbchain_mappings = defaultdict(list)
    with gzip.open(args.input, 'rt') as ipf:
        for l in ipf:
            fields = l.strip().split()
            if len(fields[0]) != 4:
                continue
            uniprot_pdbchain_mappings[fields[2]].append(fields[0] + fields[1])

    # write to a JSON file
    if args.output:
        output_file = args.output
    else:
        output_file = os.path.splitext(args.input)[0] + '.json'

    #
    with open(output_file, 'wt') as opf:
        json.dump(uniprot_pdbchain_mappings, opf, sort_keys=True, indent=4)


 
if __name__ == '__main__':
    main()
