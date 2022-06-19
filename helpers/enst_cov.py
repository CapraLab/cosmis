#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 22:45:52 2020

@author: bian
"""

import gzip
import csv
import sys
import json
from Bio import SeqIO
from argparse import ArgumentParser
from cosmis.mapping.ensembl_uniprot_pdb import SIFTS


def parse_cmd_args():
    """
    Specifies arguments that users can/should give on the command line.

    Returns
    -------
    ArgumentParser
        An object of type ArgumentParser that contains info about command-linet
    arguments.

    """
    parser = ArgumentParser(
        description='''For a given list of Ensembl transcript IDs and their
        corresponding PDB chains (identified by five-letter PDB IDs), this
        scripts computes the percentage of the residues of the transcript
        that is resolved in the PDB structure.'''
    )
    parser.add_argument(
        '-i', '--input', dest='input', required=True, type=str,
        help='''A gzip compressed file containing GENCODE annotated protein-
        coding transcript translation sequences.'''
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=True, type=str,
        help='Name of the disk file to store the computed percentages.'
    )
    parser.add_argument(
        '-m', '--mapping-list', dest='mapping_list', required=True, type=str,
        help='''A list of Ensembl transcript IDs and their corresponding 
        PDB chains (identified by five-letter PDB IDs).'''       
    )
    parser.add_argument(
        '-f', '--format', dest='format', required=False, type=str, 
        default='fasta', help='Format of the sequence file.'   
    )
    parser.add_argument(
        '-p', '--pdb-dir', dest='pdb_dir', required=True, type=str,
        help='''Path to the dictory root where PDB and SIFTS xml files are
        stored.'''
    )
    parser.add_argument(
        '-s', '--sifts-mapping', dest='sifts_mapping', required=True, 
        type=str, help='''Path to the pdb_chain_uniprot.tsv.gz file 
        downloaded from SIFTS.'''
    )
    
    return parser.parse_args()


def update_progress(completed, total):
    """
    Write progress bar to standard output.

    Parameters
    ----------
    completed : int
        The number of completed tasks.
    total : int
        Total number of tasks.

    Returns
    -------
    None.

    """
    bar_width = 50
    progress = completed / total
    n_blocks = int(round(bar_width * progress))
    text = '\rProgress: [{0}] {1}/{2}\n'.format(
        '#' * n_blocks + '-' * (bar_width - n_blocks), completed, total
    )
    sys.stdout.write(text)
    sys.stdout.flush()
    
    
def main():
    # get all command-line arguments
    cmd_args = parse_cmd_args()
    
    # read pairs of transcript to PDB chain mapping into a dictionary
    with open(cmd_args.mapping_list, 'rt') as ipf:
        enst_to_pdb = json.load(ipf)
    
    # read the transcript translation sequences into a dictionary
    translated_seqs = {}
    with gzip.open(cmd_args.input, 'rt') as handle:
        for record in SeqIO.parse(handle, cmd_args.format):
            # split record ID into fields
            id_fields = record.id.split('|')
            # the second field is the transcript ID, last field is length
            enst = id_fields[1].split('.')[0] # retain the major part
            
            try:
                seq_length = int(id_fields[-1])
            except:
                print('Invald length found in', enst, 'record!')
                continue
            
            translated_seqs[enst] = (seq_length, record.seq)

    # propress bar
    print('Now computing coverages ...')
    
    # compute the percentage of transcript residues resolved in structure
    coverages = []
    total_n_items = len(enst_to_pdb)
    completed = 0
    for k, v in enst_to_pdb.items():
        # get four-letter PDB IDs and chain IDs
        pdb_id = v[0]
        chain_id = v[1] # some PDB chain IDs could have two characters

        if not pdb_id:
            print('PDB ID is empty for {}'.format(k))
            continue
        
        # get residue-level mapping between PDB and transcript sequences
        sifts_mapping = SIFTS(cmd_args.sifts_mapping, cmd_args.pdb_dir)
        mapping = sifts_mapping.pdb_to_uniprot(pdb_id, chain_id)
        
        if mapping is None:
            print('No mapped length obtained for', pdb_id + chain_id)
            continue
        pdb_mapped_length = len(mapping)
        
        # get the length of each translated transcript sequence
        try:
            enst_length = translated_seqs[k][0]
        except KeyError:
            enst_length = None
        
        # compute and store the coverage
        try:
            coverage = pdb_mapped_length / enst_length
        except ZeroDivisionError:
            print('Zero mapped residues in', pdb_id + chain_id)
            continue
        except TypeError:
            print(
                'Failed to retrieve residue-level mapping for', 
                 pdb_id + chain_id
            )
            print('Please go to xxx to check if the file xxx is available.')
            continue
        coverages.append(
            [k, ''.join(v), enst_length, pdb_mapped_length, '%.3f' % coverage]
        )
        
        # update progress bar
        completed += 1
        update_progress(completed, total_n_items)
        
    # write the coverages to disk file201094LBzx,./
    
    print('Done computing all coverages, now writing results to file.')
    with open(cmd_args.output, 'wt') as opf:
        csv_writer = csv.writer(opf)
        csv_writer.writerow(
            [
                'transcript', 'pdb_chain', 'transcript_residues', 
                'pdb_residues', 'coverage'
            ]
        )
        for l in coverages:
            csv_writer.writerow(l)
            
    
if __name__ == '__main__':
    main()        
