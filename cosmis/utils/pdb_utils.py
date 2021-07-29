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
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import NeighborSearch
from Bio.PDB import is_aa
from cosmis.pdb_struct.contact import Contact


def get_pdb_chain(pdb_id, pdb_chain, pdb_db=None):
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
    pdb_file = os.path.join(
        pdb_dir, pdb_id + '.pdb'
    )
    if not os.path.exists(pdb_file):
        if not os.path.exists(pdb_dir):
            os.mkdir(pdb_dir)

        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(
            pdb_id, file_format='pdb',
            pdir=pdb_dir
        )

    # read in the PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    try:
        structure = pdb_parser.get_structure(id=pdb_id, file=pdb_file)
    except (FileNotFoundError, ValueError) as e:
        print('PDB file cannot be retrieved', pdb_id)
        return None
    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        print('No chain ' + pdb_chain + ' was found in ' + pdb_file)
        return None
    return chain


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

    Parameters
    ----------
    pdb_id

    Returns
    -------

    """
    pdbl = PDBList()


    if pdb_path is None or not os.path.exists(pdb_path):
        pdb_path = '/tmp/'

    # retrieve the MMCIF file for the requested PDB ID
    mmcif_file = os.path.join(pdb_path, pdb_id[1:3], pdb_id + '.cif')
    if os.path.exists(mmcif_file):
        print(mmcif_file, 'already exist. Skip downloading.')
    else:
        mmcif_file = pdbl.retrieve_pdb_file(
            pdb_id, file_format='mmCif', pdir=os.path.join(pdb_path, pdb_id[1:3])
        )

    # make a dict from meta info keys to values
    try:
        mmcif_dict = MMCIF2Dict(mmcif_file)
    except:
        return None

    # retrieve resolution
    try:
        res = mmcif_dict['_refine.ls_d_res_high']
    except KeyError:
        return None

    # remove the MMCIF file
    # os.remove(mmcif_file)

    try:
        return float(res)
    except (ValueError, TypeError):
        return None


def search_for_all_contacts(residues, radius=8):
    """
    Search for all contacts in the given set of residues based on
    distances between CB atoms.

    Parameters
    ----------
    residues
    radius

    Returns
    -------

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
    model
    atom

    Returns
    -------

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


