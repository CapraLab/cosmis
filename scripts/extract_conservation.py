#!/usr/bin/env python3

import os
import json
import csv
import sys
import gzip
import numpy as np
from argparse import ArgumentParser


def parse_cmd_args():
    """
    

    Returns
    -------
    None.

    """
    cmd_parser = ArgumentParser()
    cmd_parser.add_argument(
        '--variant-ids', '-v', dest='variant_ids', type=str, required=True,
        help='''Input IDs of variants/positions for which the conservation
        will be extracted.'''
    )
    cmd_parser.add_argument(
        '--id-type', '-t', dest='id_type', type=str, required=False,
        default='UniProt', choices=['UniProt', 'Ensembl'],
        help='Type of the IDs used in the score database file.'
    )
    cmd_parser.add_argument(
        '--database', '-d', dest='database', type=str, required=True,
        help='''Database of conservation scores keyed by Ensembl ENST IDs
        in JSON format.'''
    )
    cmd_parser.add_argument(
        '--id-mapping', '-m', dest='id_mapping', type=str, required=False,
        help='''Mapping from UniProt IDs to Ensembl ENST IDs in JSON format.'''
    )
    cmd_parser.add_argument(
        '--score', '-s', dest='score', type=str, required=False, default='ConSurf',
        choices=['GERP', 'phyloP', 'phastCons', 'RSA', 'ConSurf'],
        help='''Type of the conservation score.'''
    )
    cmd_parser.add_argument(
        '--output', '-o', dest='output', type=str, required=True,
        help='''Name of the file to which extracted scores will be written.'''
    )
    return cmd_parser.parse_args()


def main():
    """
    """
    # parse command-line arguments
    cmd_args = parse_cmd_args()
    
    # check input ID type
    if cmd_args.id_type.lower() == 'ensembl':
        if cmd_args.id_mapping is None:
            print('Please provide UniProt to Ensembl mapping file!')
            sys.exit(1)
        else:
            with open(cmd_args.id_mapping, 'rt') as ipf:
                id_mappings = json.load(ipf)
    
    # read IDs of input variants/positions
    with open(cmd_args.variant_ids, 'rt') as ipf:
        variant_ids = [l.strip() for l in ipf]
        
    # load conservation score database
    with gzip.open(cmd_args.database, 'rt') as ipf:
        conservation_db = json.load(ipf)
    
    # extract scores
    num_missing = 0
    extracted_scores = []
    for variant in variant_ids:
        acc_id, pos = variant.split('_')
        right_id = acc_id
        if cmd_args.id_type.lower() == 'ensembl':
            try:
                enst_ids = id_mappings[acc_id]
            except KeyError:
                num_missing += 1
                extracted_scores.append([acc_id, pos, 'NA'])
                continue
            # determine the right Ensembl ID
            right_id = None
            for enst_id in enst_ids:
                if enst_id in conservation_db:
                    right_id = enst_id
                    break
        # extract scores
        score_type = cmd_args.score.lower()
        if score_type in {'rsa', 'consurf'}:
            try:
                extracted_scores.append(
                    [acc_id, pos, conservation_db[right_id][int(pos) - 1]]
                )
            except (KeyError, IndexError):
                num_missing += 1
                extracted_scores.append([acc_id, pos, 'NA'])
        else:
            try:
                cons_scores = conservation_db[right_id][score_type]
                extracted_scores.append(
                    [acc_id, pos, np.mean(cons_scores[int(pos) - 1][3])]
                )
            except (KeyError, IndexError, ValueError, TypeError):
                num_missing += 1
                extracted_scores.append([acc_id, pos, 'NA'])
    
    # write extracted scores to file
    with open(os.path.abspath(cmd_args.output), 'wt') as opf:
        csv_writer = csv.writer(opf)
        csv_writer.writerow(['enst_id', 'enst_pos', 'score'])
        csv_writer.writerows(extracted_scores)
    
    print(num_missing, 'scores missing!')
    
    
if __name__ == '__main__':
    main()
