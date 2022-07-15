#!/usr/bin/env python3

import os
import sys
import json
import logging
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP

from cosmis.mapping.ensembl_uniprot_pdb import SIFTS


DSSP_BIN='/path/to/mkdssp'


def parse_cmd_args():
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed command-line arguments.
    """
    cmd_parser = ArgumentParser()
    cmd_parser.add_argument(
        '--input', '-i', dest='input', type=str, required=True,
        help='''PDB file or a file of mappings from UniProt IDs to PDB 
        files.'''
    )
    cmd_parser.add_argument(
        '--output', '-o', dest='output', type=str, required=True,
        help='''Output file containing RSAs computed by DSSP.'''
    )
    cmd_parser.add_argument(
        '--source', '-s', dest='source', type=str, required=False,
        default='PDB', choices=['PDB', 'SWISS-MODEL', 'AlphaFold'],
        help='''Source of the structure file.'''
    )
    cmd_parser.add_argument(
        '--struct-path', '-p', dest='struct_path', type=str, required=False,
        default='./pdbs/', help='''Path to the directory where all structures
        are located.'''
    )
    cmd_parser.add_argument(
        '--sifts-mapping', '-m', dest='sifts_mapping', type=str, required=False,
        default='./pdb_chain_uniprot.tsv.gz', help='''SIFTS residue-level mapping
        between UniProt sequence and PDB structures..'''
    )
    cmd_parser.add_argument(
        '--log', '-l', dest='log', type=str, required=False,
        default='uniprot_to_rsa.lg', help='''File to write logging information.'''
    )
    return cmd_parser.parse_args()


def main():
    # parse command-line arguments
    cmd_args = parse_cmd_args()

    # configure the logging system
    logging.basicConfig(
        filename=cmd_args.log,
        level=logging.INFO,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )
    
    # parse all 
    with open(cmd_args.input, 'rt') as ipf:
        uniprot_to_struct = [l.strip().split(' ') for l in ipf]
    
    uniprot_to_rsa = defaultdict(list)
    for uniprot, pep_len, struct in uniprot_to_struct:
        if cmd_args.source.lower() == 'pdb':
            struct_file = os.path.join(cmd_args.struct_path, struct[1:3], struct[:4] + '.cif')
            chain_id = struct[4:]
            cif_parser = MMCIFParser(QUIET=True)
            parsed_struct = cif_parser.get_structure('tmp', struct_file)
            
        else: 
            pdb_parser = PDBParser(PERMISSIVE=1)
            if cmd_args.source.lower() == 'swiss-model':
                struct_file = os.path.join(cmd_args.struct_path, struct + '.pdb')
                chain_id = struct[-1]
            else:
                struct_file = os.path.join(cmd_args.struct_path, struct)
                chain_id = 'A'
            parsed_struct = pdb_parser.get_structure('tmp', struct_file)
        model_for_dssp = parsed_struct[0]
        
        # call mkdssp
        print(f'Computing RSA for {uniprot} using structure at {struct_file}.')
        try:
            dssp = DSSP(model_for_dssp, struct_file, dssp=DSSP_BIN)
        except:
            logging.critical(f'Failed to compute RSA for {uniprot}.')
            continue
        
        # create a mapping of pdb_pos to rsa
        pos_to_rsa = {}
        for k, v in dssp.property_dict.items():
            if k[0] == chain_id:
                try:
                    pos_to_rsa[k[1][1]] = round(v[3], 3)
                except TypeError:
                    pos_to_rsa[k[1][1]] = np.nan

        # get residue-level mapping between PDB and transcript sequences
        if cmd_args.source.lower() == 'pdb':
            sifts_mapping = SIFTS(cmd_args.sifts_mapping, cmd_args.struct_path)
            mapping = sifts_mapping.uniprot_to_pdb(uniprot, struct[:-1], chain_id)
        else:
            mapping = {x: x for x in range(1, int(pep_len) + 1)}
        for uniprot_pos in range(1, int(pep_len) + 1):
            try:
                # map uniprot_pos to pdb_pos
                pdb_pos = mapping[uniprot_pos]
                rsa = pos_to_rsa[pdb_pos]
                uniprot_to_rsa[uniprot].append(rsa)
            except KeyError:
                uniprot_to_rsa[uniprot].append(np.nan)
                
    # write to file
    with open(cmd_args.output, 'wt') as opf:
        json.dump(uniprot_to_rsa, opf, indent=4)


if __name__ == '__main__':
    main()
