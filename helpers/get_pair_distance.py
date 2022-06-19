#!/usr/bin/env python3

from Bio.PDB import PDBParser
from argparse import ArgumentParser
import re
import csv


def parse_cmd_args():
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed command-line arguments.

    """
    parser = ArgumentParser(description='''
    ''')
    parser.add_argument(
        '-p', '--pdb', dest='pdb', required=True, type=str,
        help='A PDB file on which the distance calculation will be based.'
    )
    parser.add_argument(
        '-i', '--input', dest='input', required=True, type=str,
        help='A list of residue pairs.'
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=True, type=str,
        help='File where to write the pair distances.'
    )
    return parser.parse_args()


def get_residue(structure, chain_id, res_id):
    """
    Get the requested residue at the given chain from the given structure.

    Parameters
    ----------
    structure :
    chain_id : str
        ID of the chain.
    res_id : str
        Sequence identifier of the residue.

    Returns
    -------

    """
    chain = structure[chain_id]
    return chain[res_id]


def convert_res_id(res_id):
    """
    Convert the given residue ID to Biopython recognizable residue ID.

    Parameters
    ----------
    res_id : str, int
        A str object or an int object.

    Returns
    -------
    tuple
        A residue id represented by a tuple with three elements.

    """
    if isinstance(res_id, int):
        return ' ', res_id, ' '
    else:
        try:
            int_id = int(res_id)
            return ' ', int_id, ' '
	# non amino acid residue
        except ValueError:
            int_id, res_name = res_id.split('_')
            return 'H_' + res_name.upper(), int(int_id), ' '


def get_shortest_distance(res_a, res_b):
    """
    Compute the distances between all pairs of atoms between the given
    two residues and return the shortest.

    Parameters
    ----------
    res_a : Residue
        A Biopython Residue object.
    res_b : Residue
        A Biopython Residue object.

    Returns
    -------
    float
        The shortest distance between the given residue pair.

    """
    all_distances = []
    for a in res_a.get_atoms():
        for b in res_b.get_atoms():
            if a.get_name().startswith('H') or b.get_name().startswith('H'):
                continue
            all_distances.append(a - b)
    shortest_distance = all_distances[0]
    for dist in all_distances[1:]:
        if dist < shortest_distance:
            shortest_distance = dist
    return shortest_distance


def main():
    # parse command-line arguments
    args = parse_cmd_args()

    # parse the protein structure and get the first model
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='tmp', file=args.pdb)[0]

    # parse the residue pairs
    with open(args.input, 'rt') as ipf:
        res_pair_ids = [re.split(' ', line.strip()) for line in ipf]

    # compute pair distances
    pair_distances = []
    for a, b in res_pair_ids:
        a_id = convert_res_id(a)
        b_id = convert_res_id(b)
        res_a = get_residue(structure, 'A', a_id)
        res_b = get_residue(structure, 'A', b_id)
        pair_distances.append((a, b, get_shortest_distance(res_a, res_b)))

    # write output
    with open(args.output, 'wt') as opf:
        csv_writer = csv.writer(opf, delimiter=',')
        for pair_dist in pair_distances:
            csv_writer.writerow(pair_dist)


if __name__ == '__main__':
    main()
