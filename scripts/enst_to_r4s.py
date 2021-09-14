#!/usr/bin/env python3

import os
import json
import logging
from argparse import ArgumentParser
from collections import defaultdict


def parse_cmd_args():
    """
    Parse command-line arguments.

    Returns
    -------
    ArgumentParser
        Parse command-line arguments in an ArgumentParser object.

    """
    parser = ArgumentParser(description='''To be added ...''')
    parser.add_argument(
        '-i', '--input', dest='transcripts', type=str, required=True,
        help='''A lists of Ensembl transcript IDs for which the alignment is to
        be extracted.'''
    )
    parser.add_argument(
        '-d', '--db-path', dest='db_path', type=str, required=True,
        help='''Path to the directory where Rate4Site results are located'''
    )
    parser.add_argument(
        '-s', '--suffix', dest='suffix', type=str, required=False,
        default='_orig_rates.txt', help='''Suffix of Rate4Site result files.'''
    )
    parser.add_argument(
        '-o', '--outfile', dest='outfile', type=str, required=True,
        help='''File where to write all the parse Rate4Site results.'''
    )
    parser.add_argument(
        '--log', dest='log', type=str, required=False, default='enst_to_r4s.log',
        help='''Log file to store messages from running this script.'''
    )
    return parser.parse_args()


def main():
    # parse command-line arguments
    args = parse_cmd_args()

    # configure the logger
    logging.basicConfig(
        filename=args.log,
        level=logging.CRITICAL,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )

    # parse all the transcript ids
    with open(args.transcripts, 'rt') as ipf:
        enst_ids = [line.strip() for line in ipf]

    enst_r4s = defaultdict(list)
    print('Parsing rate4site results.')
    for enst_id in enst_ids:
        r4s_file = os.path.join(args.db_path, enst_id + args.suffix)
        if not os.path.exists(r4s_file):
            logging.critical(
                '{} does not exist.'.format(r4s_file)
            )
            continue
        with open(r4s_file, 'rt') as ipf:
            for line in ipf:
                line = line.strip()
                if line and not line.startswith('#'):
                    fields = line.split()
                    enst_r4s[enst_id].append(float(fields[2]))

    # write parsed Rate4Site results to json file
    with open(args.outfile, 'wt') as opf:
        json.dump(enst_r4s, fp=opf, indent=4)


if __name__ == '__main__':
    main()
