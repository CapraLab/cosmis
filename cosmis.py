#!/usr/bin/env python3

"""
To add a brief summary of this script.
"""

import csv
import gzip
import json
import logging
import os
import warnings
from argparse import ArgumentParser
from collections import defaultdict

import numpy as np
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.PDB import is_aa
from Bio.SeqUtils import seq1

from cosmis.mapping.sifts import SIFTS
from cosmis.utils import pdb_utils, seq_utils

warnings.simplefilter('ignore', BiopythonWarning)


def parse_cmd():
    """
    Parses command-line arguments.

    Returns
    -------
    ArgumentParser
    An object of type ArgumentParger containing parsed command-line arguments.

    """
    parser = ArgumentParser()
    parser.add_argument(
        '-c', '--config', dest='config', required=True, type=str, 
        help='''A JSON file specifying options.'''
    )
    parser.add_argument(
        '-u', '--uniprot-ids', dest='uniprot_ids', type=str, required=True, 
        help='''A list of UniProt IDs, one per line.'''
    )
    parser.add_argument(
        '-r', '--radius', dest='radius', type=float, default=8,
        help='''Radius within which to include sites.'''
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=False, type=str, 
        help='''Output file to which site-specific COSMIS values will be 
        written.'''
    )
    parser.add_argument(
        '-p', '--pdb', dest='pdb', type=str, required=False, 
        help='''PDB file representing the 3D structure of the protein.'''
    )
    parser.add_argument(
        '--multimer', dest='multimer', action='store_true', default=False,
        help='Is the input PDB file a multimer?'
    )
    parser.add_argument(
        '-w', '--overwrite', dest='overwrite', required=False, action='store_true', 
        help='''Whether to overwrite already computed COSMIS scores.'''
    )
    parser.add_argument(
        '-v', '--verbose', dest='verbose', required=False, action='store_true', 
        help='''Whether to output verbose data: number of contacting residues 
        and number of missense and synonymous variants in the neighborhood of 
        the mutation site.'''
    )
    parser.add_argument(
        '-l', '--log', dest='log', default='cosmis.log',
        help='''The file to which to write detailed computing logs.'''    
    )
                    
    return parser.parse_args()


def get_ensembl_accession(record):
    """

    Parameters
    ----------
    record : str
        Record ID is of the format: ">CCDS2.2|Hs109|chr1"

    Returns
    -------

    """
    parts = record.id.split('.')
    return parts[0]


def get_uniprot_accession(record):
    """

    Parameters
    ----------
    record

    Returns
    -------

    """
    parts = record.id.split('|')
    return parts[1]


def parse_config(config):
    """

    Parameters
    ----------
    config

    Returns
    -------

    """
    with open(config, 'rt') as ipf:
        configs = json.load(ipf)

    # do necessary sanity checks before return
    return configs


def count_variants(variants):
    """
    Collects the statistics about position-specific counts of missense and
    synonymous variants.

    Parameters
    ----------
    variants : list
        A list of variant identifiers: ['A123B', 'C456D']

    Returns
    -------
    dict
        A dictionary where the key is amino acid position and the value is
        the number of variants at this position. One dictionary for missense
        variants and one dictionary for synonymous variants.

    """
    #
    missense_counts = defaultdict(int)
    synonymous_counts = defaultdict(int)
    for variant in variants:
        vv, _, _ = variant
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        pos = vv[1:-1]  # position in the protein sequence
        # only consider rare variants
        # if int(ac) / int(an) > 0.001:
        #    continue
        if w != v:  # missense variant
            missense_counts[int(pos)] += 1
        else:  # synonymous variant
            synonymous_counts[int(pos)] += 1
    return missense_counts, synonymous_counts


def retrieve_data(uniprot_id, enst_ids, pep_dict, cds_dict, variant_dict):
    """
    """
    pep_seq = pep_dict[uniprot_id]

    valid_ensts = []
    for enst_id in enst_ids:
        cds_seq = cds_dict[enst_id].seq
        # skip if the CDS is incomplete
        if not seq_utils.is_valid_cds(cds_seq):
            print('Error: Invalid CDS.'.format(enst_id))
            continue
        if len(pep_seq) == len(cds_seq) // 3 - 1:
            valid_ensts.append(enst_id)
    if not valid_ensts:
        raise ValueError(
            'Error: {} are not compatible with {}.'.format(enst_ids, uniprot_id)
        )

    if len(valid_ensts) == 1:
        enst_id = valid_ensts[0]
        cds = cds_dict[enst_id].seq
        if enst_id not in variant_dict.keys():
           raise KeyError('Error: No record for {} in gnomAD.'.format(uniprot_id))
        variants = variant_dict[enst_id]['variants']
        return enst_id, pep_seq, cds, variants

    # if multiple transcripts are valid 
    # get the one with most variable positions
    max_len = 0
    right_enst = ''
    for enst_id in valid_ensts:
        try:
            var_pos = len(variant_dict[enst_id]['variants'])
        except KeyError:
            continue
        if max_len < var_pos:
            max_len = var_pos
            right_enst = enst_id
    cds = cds_dict[right_enst].seq
    variants = variant_dict[right_enst]['variants']

    return right_enst, pep_seq, cds, variants


def get_dataset_headers():
    """
    Returns column name for each feature of the dataset. Every time a new
    features is added, this function needs to be updated.

    Returns
    -------

    """
    header = [
        'uniprot_id', 'enst_id', 'uniprot_pos', 'uniprot_aa', 'pdb_pos',
        'pdb_aa', 'pdb_id', 'chain_id', 'seq_separations', 'num_contacts',
        'syn_var_sites', 'total_syn_sites', 'mis_var_sites', 'total_mis_sites',
        'cs_syn_poss', 'cs_mis_poss', 'cs_gc_content', 'cs_syn_prob',
        'cs_syn_obs', 'cs_mis_prob', 'cs_mis_obs', 'mis_pmt_mean', 'mis_pmt_sd',
        'mis_p_value', 'syn_pmt_mean', 'syn_pmt_sd', 'syn_p_value',
        'enst_syn_obs', 'enst_mis_obs', 'enst_syn_exp', 'enst_mis_exp', 'uniprot_length',
        'cosmis'
    ]
    return header


def main():
    # parse command-line arguments
    args = parse_cmd()

    # configure the logging system
    logging.basicConfig(
        filename=args.log,
        level=logging.INFO,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )

    # parse configuration file
    configs = parse_config(args.config)
    logging.info('Supplied configuration:')
    logging.info(json.dumps(configs, sort_keys=True, indent=4))

    # reading ENSEMBL cds
    print('Reading ENSEMBL CDS database ...')
    with gzip.open(configs['ensembl_cds'], 'rt') as cds_handle:
        cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )

    # UniProt protein sequences
    print('Reading UniProt protein sequence database ...')
    with gzip.open(configs['uniprot_pep'], 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_uniprot_accession
        )

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        variant_dict = json.load(variant_handle)

    # parse the file that maps protein UniProt IDs to PDB IDs
    with open(configs['uniprot_to_pdb'], 'rt') as ipf:
        uniprot_to_pdb = json.load(ipf)

    # parse the file that maps Ensembl transcript IDs to PDB IDs 
    with open(configs['uniprot_to_enst'], 'rt') as ipf:
        uniprot_to_enst = json.load(ipf)

    # get transcript mutation probabilities and variant counts
    print('Reading transcript mutation probabilities and variant counts ...')
    enst_mp_counts = seq_utils.read_enst_mp_count(configs['enst_mp_counts'])

    # get the directory where all output files will be stored
    output_dir = os.path.abspath(configs['output_dir'])

    # create SIFTS mapping table
    sifts_residue_mapping = SIFTS(configs['sifts_uniprot'], configs['pdb_dir'])

    # compute the COSMIS scores for each transcript
    with open(args.uniprot_ids, 'rt') as ipf:
        for line in ipf:
            uniprot_id = line.strip()
            print('Processing UniProt ID: %s' % uniprot_id)
            cosmis = []
            if args.multimer:
                output_suffix = '_cosmis_multimer.tsv'
            else:
                output_suffix = '_cosmis_monomer.tsv'
            cosmis_file = os.path.join(output_dir, uniprot_id + output_suffix)
            # skip if it was already computed and overwrite not requested
            if os.path.exists(cosmis_file) and not args.overwrite:
                print('Scores for %s already exist, skipped.' % uniprot_id)
                continue
            #
            try:
                enst_ids = uniprot_to_enst[uniprot_id]
            except KeyError:
                logging.critical(
                    'No transcript IDs were mapped to {}.'.format(uniprot_id)
                )
                continue
            try:
                right_enst, pep_seq, cds, variants = retrieve_data(
                    uniprot_id, enst_ids, pep_dict, cds_dict, variant_dict
                )
            except ValueError:
                logging.critical('No valid CDS found for {}.'.format(uniprot_id))
                continue
            except KeyError:
                logging.critical('No transcript record found for {} in gnomAD.'.format(uniprot_id))
                continue

            # calculate expected counts for each codon
            cds = cds[:-3]  # remove the stop codon
            codon_mutation_rates = seq_utils.get_codon_mutation_rates(cds)
            all_cds_ns_counts = seq_utils.count_poss_ns_variants(cds)
            all_cds_ns_sites = seq_utils.count_ns_sites(cds)

            # tabulate variants at each site
            missense_counts, synonymous_counts = count_variants(variants)

            # convert variant count to site variability
            site_variability_missense = {
                pos: 1 for pos, _ in missense_counts.items()
            }
            site_variability_synonymous = {
                pos: 1 for pos, _ in synonymous_counts.items()
            }

            # compute the total number of missense variants
            try:
                total_exp_mis_counts = enst_mp_counts[right_enst][-1]
                total_exp_syn_counts = enst_mp_counts[right_enst][-2]
            except KeyError:
                print(
                    'Transcript {} not found in {}'.format(
                        right_enst, configs['enst_mp_counts'])
                )
                continue

            # permutation test
            codon_mis_probs = [x[1] for x in codon_mutation_rates]
            codon_syn_probs = [x[0] for x in codon_mutation_rates]
            mis_p = codon_mis_probs / np.sum(codon_mis_probs)
            syn_p = codon_mis_probs / np.sum(codon_syn_probs)
            mis_pmt_matrix = seq_utils.permute_variants(
                total_exp_mis_counts, len(pep_seq), mis_p
            )
            syn_pmt_matrix = seq_utils.permute_variants(
                total_exp_syn_counts, len(pep_seq), syn_p
            )

            # get the PDB ID and PDB chain associated with this transcript
            try:
                pdb_id, pdb_chain = uniprot_to_pdb[uniprot_id]
            except KeyError:
                logging.critical(
                    '%s not found in given database %s',
                    uniprot_id, configs['uniprot_to_pdb']
                )
                logging.critical('%s was skipped ...', uniprot_id)
                print('%s not found in given database %s' %
                      (uniprot_id, configs['uniprot_to_pdb']))
                print('%s was skipped ...' % uniprot_id)
                continue

            # print message
            print(
                'Estimating COSMIS scores for:', uniprot_id, 
                right_enst, pdb_id, pdb_chain
            )

            chain = pdb_utils.get_pdb_chain(
                pdb_id, pdb_chain, configs['pdb_dir'], 'mmCif'
            )

            if chain is None:
                print('No chain was available for %s' % uniprot_id)
                continue
            
            # get all contact pairs in the PDB structure
            if args.multimer:
                structure = pdb_utils.get_structure(pdb_id, configs['pdb_dir'], 'mmCif')
                all_aa_residues = [aa for aa in structure.get_residues() if is_aa(aa, standard=True)]
            else:
                all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa, standard=True)]
            all_contacts = pdb_utils.search_for_all_contacts(
                all_aa_residues, radius=args.radius
            )

            # index all contacts by residue ID
            indexed_contacts = defaultdict(list)
            for c in all_contacts:
                indexed_contacts[c.get_res_a()].append(c.get_res_b())
                indexed_contacts[c.get_res_b()].append(c.get_res_a())

            # one-to-one mapping between UniProt residues and PDB residues
            # used later in determining count of variants in contact set
            uniprot_to_pdb_mapping = sifts_residue_mapping.uniprot_to_pdb(
                uniprot_id, pdb_id, pdb_chain
            )
            pdb_to_uniprot_mapping = {v: k for k, v in uniprot_to_pdb_mapping.items()}
            if uniprot_to_pdb_mapping is None or pdb_to_uniprot_mapping is None:
                print('Mapping between UniProt and PDB failed! Skipped!')
                continue

            # row identifier for each COSMIS score record
            id_fields = [uniprot_id, right_enst]
            
            for i, a in enumerate(pep_seq, start=1):
                # @TODO need to get the corresponding position in the PDB file
                # It is highly likely that the amino acid at position i in the
                # PDB file will not be the same amino acid at position i in the
                # ENSP peptide sequence
                try:
                    pdb_pos = uniprot_to_pdb_mapping[i]
                    res = chain[pdb_pos]
                except KeyError:
                    logging.critical(
                        'Residue %s in %s not found in chain %s in PDB file: %s', 
                        i, uniprot_id, pdb_chain, pdb_id
                    )
                    continue
                    
                # check that the amino acid in ENSP sequence matches 
                # that in the PDB structure
                pdb_aa = seq1(res.get_resname())
                if a != pdb_aa:
                    logging.critical(
                        'Residue in %s did not match that in PDB %s '
                        'chain %s at %s: %s vs %s',
                        uniprot_id, pdb_id, pdb_chain, i, a, pdb_aa
                    )
                    continue

                contact_res = indexed_contacts[res]
                # only consider contacts within the same protein
                contacts_pdb_pos = [
                    r.get_id()[1] for r in contact_res if r.get_id()[1] in pdb_to_uniprot_mapping.keys()
                ]
                seq_seps = ';'.join(str(x) for x in
                                    [pdb_to_uniprot_mapping[c] - i for c in contacts_pdb_pos])
                num_contacts = len(contacts_pdb_pos)

                total_missense_obs = missense_counts.setdefault(i, 0)
                total_synonymous_obs = synonymous_counts.setdefault(i, 0)
                mis_var_sites = site_variability_missense.setdefault(i, 0)
                syn_var_sites = site_variability_synonymous.setdefault(i, 0)
                total_missense_poss = all_cds_ns_counts[i - 1][0]
                total_synonyms_poss = all_cds_ns_counts[i - 1][1]
                total_missense_sites = all_cds_ns_sites[i - 1][0]
                total_synonymous_sites = all_cds_ns_sites[i - 1][1]
                total_synonymous_rate = codon_mutation_rates[i - 1][0]
                total_missense_rate = codon_mutation_rates[i - 1][1]

                # get sequence context of the contact set
                try:
                    seq_context = seq_utils.get_codon_seq_context(
                        [pdb_to_uniprot_mapping[x] for x in contacts_pdb_pos] + [i], cds
                    )
                except IndexError:
                    logging.critical(
                        'Error in retrieving sequence context:', 
                        contacts_pdb_pos, right_enst, uniprot_id
                    )
                    break

                # compute the GC content of the sequence context
                if len(seq_context) == 0:
                    print('No nucleotides were found in sequence context!')
                    continue
                gc_fraction = seq_utils.gc_content(seq_context)

                if contacts_pdb_pos:
                    all_uniprot_pos = []
                    for j in contacts_pdb_pos:
                        # @TODO need to get the corresponding position in ENSP
                        # amino acid sequence using the SIFTS mapping table
                        try:
                            uniprot_pos = pdb_to_uniprot_mapping[j]
                        except KeyError:
                            logging.critical(
                                'PDB residue %s in %s chain %s not found '
                                'in %s, skipped',
                                j, pdb_id, pdb_chain, uniprot_id
                            )
                            continue

                        # count the total # observed variants in contacting residues
                        total_missense_obs += missense_counts.setdefault(uniprot_pos, 0)
                        total_synonymous_obs += synonymous_counts.setdefault(uniprot_pos, 0)
                        mis_var_sites += site_variability_missense.setdefault(uniprot_pos, 0)
                        syn_var_sites += site_variability_synonymous.setdefault(uniprot_pos, 0)
                        # count the total # expected variants
                        try:
                            total_missense_poss += all_cds_ns_counts[uniprot_pos - 1][0]
                            total_synonyms_poss += all_cds_ns_counts[uniprot_pos - 1][1]
                            total_missense_sites += all_cds_ns_sites[uniprot_pos - 1][0]
                            total_synonymous_sites += all_cds_ns_sites[uniprot_pos - 1][1]
                            total_synonymous_rate += codon_mutation_rates[uniprot_pos - 1][0]
                            total_missense_rate += codon_mutation_rates[uniprot_pos - 1][1]
                        except IndexError:
                            logging.critical(
                                'PDB residue %s in %s chain %s out of the range '
                                'of %s codons, something must be wrong with SIFTS mapping',
                                j, pdb_id, pdb_chain, right_enst
                            )
                            continue
                        all_uniprot_pos.append(uniprot_pos)

                    # get permutation statistics
                    mis_pmt_mean, mis_pmt_sd, mis_p_value = seq_utils.get_permutation_stats(
                        mis_pmt_matrix, all_uniprot_pos + [i], total_missense_obs
                    )
                    syn_pmt_mean, syn_pmt_sd, syn_p_value = seq_utils.get_permutation_stats(
                        syn_pmt_matrix, all_uniprot_pos + [i], total_synonymous_obs
                    )

                    # push results for the current residue
                    cosmis.append(
                        id_fields + 
                        [i, a, pdb_pos, pdb_aa] + 
                        [pdb_id, pdb_chain] + 
                        [
                            seq_seps, 
                            num_contacts + 1,
                            syn_var_sites, 
                            '{:.3f}'.format(total_synonymous_sites),
                            mis_var_sites, 
                            '{:.3f}'.format(total_missense_sites),
                            total_synonyms_poss,
                            total_missense_poss,
                            '{:.3f}'.format(gc_fraction),
                            '{:.3e}'.format(total_synonymous_rate),
                            total_synonymous_obs, 
                            '{:.3e}'.format(total_missense_rate),
                            total_missense_obs,
                            '{:.3f}'.format(mis_pmt_mean),
                            '{:.3f}'.format(mis_pmt_sd),
                            '{:.3e}'.format(mis_p_value),
                            '{:.3f}'.format(syn_pmt_mean),
                            '{:.3f}'.format(syn_pmt_sd),
                            '{:.3e}'.format(syn_p_value),
                            enst_mp_counts[right_enst][2],
                            enst_mp_counts[right_enst][4],
                            total_exp_syn_counts,
                            total_exp_mis_counts,
                            len(pep_seq),
                            '{:.3f}'.format((total_missense_obs - mis_pmt_mean) / mis_pmt_sd)
                        ]
                    )
                else:
                    print("No contacts found.")
                    continue

            with open(file=cosmis_file, mode='wt') as opf:
                header = get_dataset_headers()
                csv_writer = csv.writer(opf, delimiter='\t')
                csv_writer.writerow(header)
                csv_writer.writerows(cosmis)


if __name__ == '__main__':
    main()
