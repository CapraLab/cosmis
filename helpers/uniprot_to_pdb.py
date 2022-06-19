#!/usr/bin/env python3

"""
    @summary For each given Ensembl transcript ID, determine the PDB ID of 
    the structure that best represents the transcript.
    @author Bian Li
    @contact bian.li@vanderbilt.edu
    @created: 08/03/2021
    @change 08/03/2021

"""

import json
from argparse import ArgumentParser
from cosmis.mapping.ensembl_uniprot_pdb import EnsemblUniProtPDB


def parse_cmd():
    """
    Parse command-line arguments
    
    Returns
    -------
    ArgumentParser
    An object of type ArgumentParser containing parsed commond-line arguments.

    """
    parser = ArgumentParser()
    parser.add_argument(
        '-u', '--uniprot_ids', dest='uniprot_ids', required=True, type=str, 
        help='A list of UniProt IDs, one per line.'
    )
    parser.add_argument(
        '--sifts-ensembl', dest='sifts_ensembl', required=True, type=str, 
        help='SIFTS chain-level mapping between UniProt, PDB, and Ensembl.'
    )
    parser.add_argument(
        '--sifts-uniprot', dest='sifts_uniprot', required=True, type=str, 
        help='SIFTS residue-level mapping between UniProt and PDB.'
    )
    parser.add_argument(
        '-d', '--database', dest='pdb_path', required=False, type=str, 
        help='Path to local PDB database.'
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=False, type=str,
        help='Output file.'
    )
    parser.add_argument(
        '--multimeric-state', dest='multimeric_state', type=int, default=0,
        help='Consider only homo-oligomers.'
    )
    args = parser.parse_args()
    
    # do any necessary check on command-line arguments here
    return args


def main():
    # parse command-line arguments
    args = parse_cmd()

    # read in the mapping file
    ensembl_pdbchain_mappings = EnsemblUniProtPDB(
        args.sifts_ensembl, args.sifts_uniprot, args.pdb_path
    )

    # get supplied Ensembl transcript IDs
    with open(args.uniprot_ids, 'rt') as ipf:
        uniprot_ids = [l.strip() for l in ipf]

    #
    uniprot_uniq_pdbchains = {}
    for uniprot_id in uniprot_ids:
        print('\nRetrieving best PDB chain for %s' % uniprot_id)
        pdb_id, pdb_chain = ensembl_pdbchain_mappings.uniprot_to_pdb(
            uniprot_id, args.multimeric_state
        )
        if pdb_id is None or pdb_chain is None:
            print('No PDB chains were found for %s' % uniprot_id)
            continue
        uniprot_uniq_pdbchains[uniprot_id] = [pdb_id, pdb_chain]

    with open(args.output, 'wt') as opf:
        json.dump(uniprot_uniq_pdbchains, opf, sort_keys=True, indent=4)


if __name__ == '__main__':
    main()
