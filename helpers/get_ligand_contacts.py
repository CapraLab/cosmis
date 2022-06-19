#!/usr/bin/env python3

from Bio.PDB import PDBParser
from argparse import ArgumentParser
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
        help='''A PDB file on which the distance calculation will be based.'''
    )
    parser.add_argument(
        '-l', '--ligand_id', dest='ligand_id', required=True, type=str,
        help='''Ligand residue id in the format: 000_XXX where XXX is the 
        three-letter ligand name, 000 is the ligand residue number.'''
    )
    parser.add_argument(
        '-c', '--chain-id', dest='chain_id', required=False, type=str,
        default='A', help='''Chain ID of the protein and ligand.'''           
    )
    parser.add_argument(
        '-t', '--threshold', dest='threshold', required=False, type=float,
        default=5.0, help='''Distance threshold.'''           
    )  
    parser.add_argument(
        '-o', '--output', dest='output', required=True, type=str,
        help='''File where to write the pair distances.'''
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

    # get the ligand
    ligand_res_id = convert_res_id(args.ligand_id)
    protein_chain_id = args.chain_id[0]
    ligand_chain_id = args.chain_id[1]
    ligand_res = get_residue(structure, ligand_chain_id, ligand_res_id)

    # compute pair distances
    pair_distances = []
    for res in structure[protein_chain_id].get_residues():
        if res == ligand_res:
            continue
        res_id = res.get_id()[1]
        shortest_distance = get_shortest_distance(ligand_res, res)
        if shortest_distance <= args.threshold:
            pair_distances.append((res_id, '%.3f' % shortest_distance))

    # write output
    with open(args.output, 'wt') as opf:
        csv_writer = csv.writer(opf, delimiter=',')
        for pair_dist in pair_distances:
            csv_writer.writerow(pair_dist)


if __name__ == '__main__':
    main()
