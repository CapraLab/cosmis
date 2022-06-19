#!/usr/bin/env python3

import json
from argparse import ArgumentParser


def parse_cmd():
    """
    Parse command-line arguments.

    Returns
    -------
    An ArgumentParser object.

    """
    # add flags and specify arguments
    cmd_parser = ArgumentParser(
        description='''To be added ...'''
    )
    cmd_parser.add_argument(
        '--input', '-i', dest='input', type=str, required=True,
        help='''A file that contains one pair of UniProt IDs on each line, representing a protein
        complex.'''
    )
    cmd_parser.add_argument(
        '--mapping', '-m', dest='mapping', type=str, required=True,
        help='''A JSON file that contains the mapping of UniProt IDs to "best" structures. Protein
        structures are represented by PDB IDs + chain IDs.'''
    )
    cmd_parser.add_argument(
        '--output', '-o', dest='output', type=str, required=True,
        help='File to store results.'
    )
    # parse and return flags and arguments
    return cmd_parser.parse_args()


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    cmd_args = parse_cmd()

    #
    with open(cmd_args.input, 'rt') as ipf:
        uniprot_pairs = [line.strip().split(':') for line in ipf]

    # load the mapping file
    with open(cmd_args.mapping, 'rt') as ipf:
        uniprot_pdb_mapping = json.load(ipf)

    # get the pdb file for each of the protein pairs
    complex_mappings = []
    for p1, p2 in uniprot_pairs:
        try:
            pdb1 = uniprot_pdb_mapping[p1]
        except KeyError:
            print(f'No PDB ID found for {p1}')
            continue
        try:
            pdb2 = uniprot_pdb_mapping[p2]
        except KeyError:
            print(f'No PDB ID found for {p2}')
            continue

        if pdb1[0] != pdb2[0]:
            print(f'{p1} and {p2} are not mapped to the same PDB file.')
            print(f'{p1}: {pdb1[0]}\t{p2}: {pdb2[0]}')
            continue

        # add the PDB chain to UniProt ID mapping for the current protein complex
        complex_mappings.append(
            f'{pdb1[0]}\t{pdb1[1]}:{p1}\t{pdb2[1]}:{p2}'
        )

    # write the mappings to file
    with open(cmd_args.output, 'wt') as opf:
        opf.write('\n'.join(complex_mappings))


if __name__ == '__main__':
    main()