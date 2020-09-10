#!/usr/bin/env python3

"""
    @summary
    @author
    @contact
    @change

"""

import os
import json
import gzip
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.Alphabet import ProteinAlphabet
from Bio.PDB import PDBParser, PPBuilder


def parse_cmd():
    """
    
    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument('-u', '--uniprot', dest='uniprot_ids', required=True,
                        type=str, help='A list of UniProt IDs, one per line.')
    parser.add_argument('-m', '--mappings', dest='sifts_mappings', required=True,
                        type=str, help='SIFTS PDB chain to UniProt mapping file.')
    parser.add_argument('-d', '--database', dest='pdb_path', required=False,
                        type=str, help='Path to local PDB database.')
    parser.add_argument('-o', '--output', dest='output', required=False, type=str,
                        help='Output file.')
    args = parser.parse_args()
    # do any necessary check on command-line arguments here
    return args


def get_pdb_chain(five_letter_id=None, pdb_path=None):
    """
    Creates a Bio.PDB.Chain object for the requested PDB chain.

    Parameters
    ----------
    five_letter_id : str
        Identifier of the PDB chain as a five-letter string.
    pdb_path : str
        Path to the local PDB database.

    Returns
    -------
    Bio.PDB.Chain
        The requested PDB chain as a Bio.PDB.Chain object.

    """
    pdb_parser = PDBParser(PERMISSIVE=1)
    pdb_id = five_letter_id[:-1]
    chain_id = five_letter_id[-1].upper()
    pdb_file = 'pdb' + pdb_id + '.ent.gz'
    path_to_pdb_file = os.path.join(pdb_path, pdb_file)
    # download pdb from the web if not already exists
    chain = None
    if os.path.exists(path_to_pdb_file):
        with gzip.open(path_to_pdb_file, 'rt') as opf:
            structure = pdb_parser.get_structure(five_letter_id[:-1], opf)
        try:
            chain = structure[0][chain_id]
        except KeyError:
            print('No chain ' + chain_id + ' was found in ' + path_to_pdb_file)
    return chain


def get_chain_seq(chain=None):
    """

    Parameters
    ----------
    chain

    Returns
    -------
    Bio.Seq

    """
    # build a polypeptide object from given chain
    ppb = PPBuilder()
    chain_seq = Seq('', ProteinAlphabet())
    for pp in ppb.build_peptides(chain):
        chain_seq += pp.get_sequence()
    return chain_seq


def get_resolution(pdb_id=None):
    """

    Parameters
    ----------
    pdb_id

    Returns
    -------

    """
    pass


def main():
    # parse command-line arguments
    args = parse_cmd()

    # read in the mapping file
    with open(args.sifts_mappings, 'rt') as ipf:
        uniprot_pdbchain_mappings = json.load(ipf)

    # get supplied UniProt IDs
    with open(args.uniprot_ids, 'rt') as ipf:
        uniprot_ids = [l.strip() for l in ipf]

    #
    with open(args.output, 'wt') as opf:
        for u in uniprot_ids:
            try:
                pdbchains = uniprot_pdbchain_mappings[u]
            except KeyError:
                print('No PDB chain found for UniProt ID: ' + u)
            unique_pdbchain = ""
            max_len = 0
            for c in pdbchains:
                chain = get_pdb_chain(c, args.pdb_path)
                if chain is not None:
                    chain_seq = get_chain_seq(chain)
                    if len(chain_seq) > max_len:
                        max_len = len(chain_seq)
                        unique_pdbchain = c
            opf.write(u + '\t' + unique_pdbchain + '\n')


if __name__ == '__main__':
    main()
