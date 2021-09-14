#!/usr/bin/env python3

"""
    @summary For each given Ensembl transcript ID, determine the PDB ID of 
    the structure that best represents the transcript.
    @author Bian Li
    @contact bian.li@vanderbilt.edu
    @change 09/07/2020

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
        '-t', '--enst_ids', dest='enst_ids', required=True, type=str, 
        help='A list of ENST IDs, one per line.'
    )
    parser.add_argument(
        '-m', '--mappings', dest='sifts_mappings', required=True, type=str, 
        help='SIFTS PDB chain to UniProt mapping file.'
    )
    parser.add_argument(
        '-d', '--database', dest='pdb_path', required=False, type=str, 
        help='Path to local PDB database.'
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=False, type=str,
        help='Output file.'
    )
    args = parser.parse_args()
    
    # do any necessary check on command-line arguments here
    return args


def main():
    # parse command-line arguments
    args = parse_cmd()

    # read in the mapping file
    ensembl_pdbchain_mappings = EnsemblUniProtPDB(
        args.sifts_mappings, path_path=args.pdb_path
    )

    # get supplied Ensembl transcript IDs
    with open(args.enst_ids, 'rt') as ipf:
        enst_ids = [l.strip() for l in ipf]

    #
    enst_uniq_pdbchains = {}
    for enst in enst_ids:
        print('Retrieving best PDB chain for %s' % enst)
        pdb_id, pdb_chain = ensembl_pdbchain_mappings.enst_to_pdb(enst)
        if pdb_id is None or pdb_chain is None:
            print('No PDB chains were found for %s' % enst)
            continue
        enst_uniq_pdbchains[enst] = [pdb_id, pdb_chain]

    with open(args.output, 'wt') as opf:
        json.dump(enst_uniq_pdbchains, opf, sort_keys=True, indent=4)


if __name__ == '__main__':
    main()
