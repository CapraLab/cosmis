#!/usr/bin/env python3

"""
    @summary
    @author
    @contact
    @change

"""

import os
import numpy as np
from Bio.Seq import Seq
from Bio.PDB import PDBParser, PDBList, PPBuilder
from Bio.PDB import parse_pdb_header
from Bio.PDB import NeighborSearch
from Bio.PDB import is_aa
from Bio.PDB import MMCIFParser
from cosmis.pdb_struct.contact import Contact


def get_pdb_chain(pdb_id, pdb_chain, pdb_db=None, pdb_format='pdb'):
    """
    Creates a Bio.PDB.Chain object for the requested PDB chain.

    Parameters
    ----------
    pdb_id : str
        Identifier of the PDB chain as a five-letter string.
    pdb_chain : str
        Identifier of the PDB chain as a five-letter string.
    pdb_db : str
        Path to the local PDB database.

    Returns
    -------
    Bio.PDB.Chain
        The requested PDB chain as a Bio.PDB.Chain object.

    """
    # check to see if the PDB file is already downloaded
    # if not, download and store it locally
    pdb_dir = os.path.join(pdb_db, pdb_id[1:3])
    if pdb_format.lower() == 'pdb':
        file_suffix = '.pdb'
    else:
        file_suffix = '.cif'
    pdb_file = os.path.join(
        pdb_dir, pdb_id + file_suffix
    )
    if not os.path.exists(pdb_file):
        if not os.path.exists(pdb_dir):
            os.mkdir(pdb_dir)

        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(
            pdb_id, file_format=pdb_format,
            pdir=pdb_dir
        )

    # read in the PDB file
    if pdb_format.lower() == 'pdb':
        pdb_parser = PDBParser(PERMISSIVE=1)
    else: 
        pdb_parser = MMCIFParser(QUIET=True)
    try:
        structure = pdb_parser.get_structure(pdb_id, pdb_file)
    except (FileNotFoundError, ValueError) as e:
        print('PDB file cannot be retrieved', pdb_id)
        return None
    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        print('No chain ' + pdb_chain + ' was found in ' + pdb_file)
        return None
    return chain


def get_structure(pdb_id, pdb_db=None, pdb_format='pdb'):
    """
    Create a Bio.PDB.Structure object for the requested PDB ID.

    Parameters
    ----------
    pdb_id : str
        Four-letter identifier of the PDB file.
    pdb_db : str
        Path to the local PDB database.
    pdb_format : str
        PDB or MMCif

    Returns
    -------
    Structure

    """
    # check to see if the PDB file is already downloaded
    # if not, download and store it locally
    pdb_dir = os.path.join(pdb_db, pdb_id[1:3])
    if pdb_format.lower() == 'pdb':
        file_suffix = '.pdb'
    else:
        file_suffix = '.cif'
    pdb_file = os.path.join(
        pdb_dir, pdb_id + file_suffix
    )
    if not os.path.exists(pdb_file):
        if not os.path.exists(pdb_dir):
            os.mkdir(pdb_dir)

        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(
            pdb_id, file_format=pdb_format,
            pdir=pdb_dir
        )

    # read in the PDB file
    if pdb_format.lower() == 'pdb':
        pdb_parser = PDBParser(PERMISSIVE=1)
    else:
        pdb_parser = MMCIFParser(QUIET=True)
    try:
        return pdb_parser.get_structure(pdb_id, pdb_file)
    except (FileNotFoundError, ValueError) as e:
        print('PDB file cannot be retrieved', pdb_id)
        return None


def get_chain_seq(chain):
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


def get_resolution(pdb_id, pdb_path=None):
    """
    Get the resolution of the structure represented by the given PDB ID.

    Parameters
    ----------
    pdb_id : str
        Four-letter identifier of the PDB file.
    pdb_path : str
        Path to the local database where the PDB files are stored.

    Returns
    -------
    float
        Resolution of the structure in angstroms.

    """
    pdbl = PDBList()

    if pdb_path is None:
        pdb_path = '/tmp/'
    if not os.path.exists(pdb_path):
        os.mkdir(pdb_path)

    # retrieve the PDB file for the requested PDB ID
    cif_file = os.path.join(pdb_path, pdb_id[1:3], pdb_id + '.cif')
    if not os.path.exists(cif_file) or os.stat(cif_file).st_size == 0:
        cif_file = pdbl.retrieve_pdb_file(
            pdb_id, file_format='mmCif', pdir=os.path.join(pdb_path, pdb_id[1:3])
        )

    # make a dict from meta info keys to values
    mmcif_parser = MMCIFParser(QUIET=True)
    structure_cif = mmcif_parser.get_structure(pdb_id, cif_file)
    try:
        res = structure_cif.header['resolution']
    except KeyError:
        return None

    try:
        return float(res)
    except (ValueError, TypeError):
        return None


def search_for_all_contacts(residues, radius=8.0):
    """
    Search for all contacts in the given set of residues based on
    distances between CB atoms.

    Parameters
    ----------
    residues : list
        A list of Biopython Residue objects.
    radius : float
        The radius within which two residues are considered in contact.

    Returns
    -------
    list
        A list of Contact objects.

    """
    atom_list = []
    for r in residues:
        if r.get_resname() == 'GLY':
            try:
                atom_list.append(r['CA'])
            except KeyError:
                print('No CA atom found for GLY:', r, 'skipped ...')
                continue
        else:
            try:
                atom_list.append(r['CB'])
            except KeyError:
                print('No CB atom found for:', r.get_resname(), 'skipped ...')
                continue
            # atom_list += [a for a in r.get_atoms() if a.get_name()
            #               not in BACKBONE_ATOMS]
        # atom_list += [a for a in r.get_atoms()]
    ns = NeighborSearch(atom_list)
    all_contacts = [
        Contact(res_a=c[0], res_b=c[1])
        for c in ns.search_all(radius, level='R')
    ]
    return all_contacts


def compute_adjacency_matrix(model, cutoff):
    """

    Parameters
    ----------
    model : Bio.PDB.Model
        A protein structure model stored as a Bio.PDB.Model object.
    cutoff : float
        Distance cutoff used to consider that two residues are adjacent.

    Returns
    -------
    NumPy 2DArray

    """
    pass


def compute_distance_matrix(model, atom='CB'):
    """
    Computes an L * L matrix consisting the distances between residues.

    Parameters
    ----------
    model : Bio.PDB.Chain
        Structure of the protein for which to compute a distance matrix.
    atom : str
        Name of the atom based on which to measure inter-residue distance.

    Returns
    -------
    NumPy ndarray
        A matrix of inter-residue distances.

    """
    # get all amino acid residues
    all_aa_res = []
    for res in model.get_residues():
        if is_aa(res):
            all_aa_res.append(res)

    # get all measurement points
    measurement_points = []
    for res in all_aa_res:
        # get the right atom
        if atom == 'CB' and res.get_resname() == 'GLY':
            point = 'CA'
        else:
            point = atom
        try:
            res_atom = res[point]
        except KeyError:
            print('Missing {} in {}. Skipped!'.format(point, res))
            continue
        measurement_points.append(res_atom)

    # compute distances
    res_distances = []
    for point_a in measurement_points:
        res_a_distances = []
        for point_b in measurement_points:
            # compute the distance between residues A and B
            a_b_distance = point_a - point_b
            res_a_distances.append(a_b_distance)
        res_distances.append(res_a_distances)
    return np.array(res_distances)


def compute_distance_dict(model, atom='CB'):
    """
    Computes an L * L matrix consisting the distances between residues.

    Parameters
    ----------
    model : Bio.PDB.Chain
        Structure of the protein for which to compute a distance matrix.
    atom : str
        Name of the atom based on which to measure inter-residue distance.

    Returns
    -------
    dict
        A dictionary of inter-residue distances, keyed by residue position pairs.

    """
    # get all amino acid residues
    all_aa_res = []
    for res in model.get_residues():
        if is_aa(res):
            all_aa_res.append(res)

    # get all measurement points
    distances = {}
    for res_a in all_aa_res:
        # get the right atom
        point_a = atom
        if atom == 'CB' and res_a.get_resname() == 'GLY':
            point_a = 'CA'
        try:
            res_a_atom = res_a[point_a]
        except KeyError:
            print('Missing {} in {}. Skipped!'.format(point_a, res_a))
            continue
        for res_b in all_aa_res:
            # get the right atom
            point_b = atom
            if atom == 'CB' and res_b.get_resname() == 'GLY':
                point_b = 'CA'
            try:
                res_b_atom = res_b[point_b]
            except KeyError:
                print('Missing {} in {}. Skipped!'.format(point_b, res_b))
                continue

            distances[(res_a.get_id()[1], res_b.get_id()[1])] = res_a_atom - res_b_atom

    return distances


def compute_mean_distance(residues):
    """

    Parameters
    ----------
    residues

    Returns
    -------

    """
    # get the atoms for computing the distances
    atoms = []
    for r in residues:
        if r.get_resname() == 'GLY':
            try:
                atoms.append(r['CA'])
            except KeyError:
                print(f'No CA atom found for GLY: {r} skipped ...')
                continue
        else:
            try:
                atoms.append(r['CB'])
            except KeyError:
                print(f'No CB atom found for: {r} skipped ...')
                continue

    # compute the mean pairwise distance
    distances = []
    for i, a1 in enumerate(atoms):
        for _, a2 in enumerate(atoms, i + 1):
            distances.append(a1 - a2)

    return np.mean(distances)
