#!/usr/bin/env python3

"""
    @summary
    @author
    @contact
    @change

"""

import os
from Bio.Seq import Seq
from Bio.Alphabet import ProteinAlphabet
from Bio.PDB import PDBParser, PDBList, PPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


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

