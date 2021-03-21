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
        '-t', '--transcripts', dest='transcripts', type=str, required=True, 
        help='''A list of ENSEMBL transcript IDs, one per line.'''
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


def get_ccds_accession(record):
    """

    Parameters
    ----------
    record

    Returns
    -------

    """
    parts = record.id.split('|')
    return parts[0]


def get_transcript_pep_seq(enst_id, ensp_id, pep_dict):
    """

    Parameters
    ----------
    enst_id : str

    ensp_id : str

    pep_dict : dict

    Returns
    -------

    """
    # if len(ensp_id) > 1:
    #     logging.critical(
    #         'More than one ENSP IDs found for transcript %s: %s',
    #         enst_id, ensp_id
    #     )
    #     sys.exit(1)
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        logging.critical('%s not found in given database', ensp_id)
        logging.critical('%s was skipped ...', enst_id)
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


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


def get_dataset_headers():
    """
    Returns column name for each feature of the dataset. Every time a new
    features is added, this function needs to be updated.

    Returns
    -------

    """
    header = [
        'enst_id', 'ensp_id', 'uniprot_id', 'ensp_pos', 'ensp_aa', 'pdb_pos',
        'pdb_aa', 'pdb_id', 'chain_id', 'seq_separations', 'num_contacts',
        'syn_var_sites', 'total_syn_sites', 'mis_var_sites', 'total_mis_sites',
        'cs_syn_poss', 'cs_mis_poss', 'cs_gc_content', 'cs_syn_prob',
        'cs_syn_obs', 'cs_mis_prob', 'cs_mis_obs', 'pmt_mean', 'pmt_sd',
        'p_value', 'phylop', 'gerp', 'r4s', 'enst_syn_obs', 'enst_mis_obs',
        'enst_mis_exp', 'ensp_length'
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
        ensembl_cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )

    # CCDS concensus coding sequences
    print('Reading NCBI CCDS database ...')
    with gzip.open(configs['ccds_cds'], 'rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(
            SeqIO.parse(ccds_handle, format='fasta'),
            key_function=get_ccds_accession
        )

    # ENSEMBL peptide sequences
    print('Reading ENSEMBL protein sequence database ...')
    with gzip.open(configs['ensembl_pep'], 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_ensembl_accession
        )

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    # parse the file that maps Ensembl transcript IDs to PDB IDs 
    with open(configs['enst_to_pdb'], 'rt') as ipf:
        enst_to_pdb = json.load(ipf)

    # get phylop scores
    print('Loading PhyloP scores ...')
    with gzip.open(configs['enst_to_phylop'], 'rt') as ipf:
        enst_to_phylop = json.load(ipf)

    # load GERP++ scores
    print('Loading GERP++ scores ...')
    with gzip.open(configs['enst_to_gerp'], 'rt') as ipf:
        enst_to_gerp = json.load(ipf)

    # load GERP++ scores
    print('Loading rate4site scores ...')
    with gzip.open(configs['enst_to_r4s'], 'rt') as ipf:
        enst_to_r4s = json.load(ipf)

    # get transcript mutation probabilities and variant counts
    print('Reading transcript mutation probabilities and variant counts ...')
    enst_mp_counts = seq_utils.read_enst_mp_count(configs['enst_mp_counts'])

    # get the directory where all output files will be stored
    output_dir = os.path.abspath(configs['output_dir'])

    # create SIFTS mapping table
    sifts_residue_mapping = SIFTS(configs['sifts_residue_mapping_file'], configs['pdb_dir'])

    # compute the COSMIS scores for each transcript
    with open(args.transcripts, 'rt') as ipf:
        for transcript in ipf:
            transcript = transcript.strip()
            print('Processing transcript %s' % transcript)
            cosmis = []
            cosmis_file = os.path.join(output_dir, transcript + '_cosmis.tsv')
            # skip if it was already computed and overwrite not requested
            if os.path.exists(cosmis_file) and not args.overwrite:
                print('Scores for %s already exist, skipped.' % transcript)
                continue
            
            # get the amino acid sequence of the transcript
            try:
                ensp_id = transcript_variants[transcript]['ensp'][0]
            except KeyError:
                logging.critical(
                    'Transcript %s not found in %s',
                    transcript, configs['gnomad_variants']
                )
                continue

            transcript_pep = get_transcript_pep_seq(
                transcript, ensp_id, pep_dict
            )

            # uniprot id for retrieving the corresponding PDB chain
            uniprot_ids = transcript_variants[transcript]['swissprot']
            if len(uniprot_ids) > 1:
                logging.critical(
                    'ERROR: more than one UniProt IDs found for transcript ' 
                    '%s: %s', transcript, ','.join(uniprot_ids)
                )
                continue
            
            # get the UniProt ID
            uniprot_id = uniprot_ids[0]

            if transcript_pep is None:
                logging.critical(
                    'No peptide sequence found for %s in %s' % 
                    (transcript, configs['ensembl_pep']) 
                )

                # fall back to UniProt sequence
                uniprot_pep = seq_utils.get_uniprot_aa_seq(uniprot_id)
                if uniprot_pep is None:
                    continue
                else:
                    transcript_pep = uniprot_pep

            # get all variants of this transcript reported in gnomAD
            try:
                variants = transcript_variants[transcript]['variants']
            except KeyError:
                logging.critical(
                    'No variants found in %s in gnomAD', transcript
                )
                logging.critical('%s was skipped ...', transcript)
                continue

            # get the coding sequence of the transcript
            try:
                transcript_cds = ensembl_cds_dict[transcript].seq
            except KeyError:
                logging.critical(
                    'No CDS found in Ensembl CDS database! Looking for it '
                    'in the CCDS database ...'                    
                )
                transcript_cds = None

            if transcript_cds is None:
                try:
                    ccds_id = transcript_variants[transcript]['ccds'][0]
                    transcript_cds = ccds_dict[ccds_id].seq
                except KeyError:
                    logging.critical(
                        'No CDS found in CCDS database! Skipped.'
                    )
                    continue

            # skip if the CDS is incomplete
            if not seq_utils.is_valid_cds(transcript_cds):
                print('Error: Invalid CDS. Shipped {}.'.format(transcript))
                continue

            # get the phyloP scores for the current transcript
            try:
                transcript_phylop_scores = enst_to_phylop[transcript]['phylop']
            except KeyError:
                print('No phyloP scores are available for {}'.format(transcript))
                continue

            # get the GERP++ scores for the current transcript
            try:
                transcript_gerp_scores = enst_to_gerp[transcript]['gerp']
            except KeyError:
                print(
                    'No GERP++ scores are available for {}'.format(transcript))
                continue

            # get the rate4site scores for the current transcript
            try:
                transcript_r4s_scores = enst_to_r4s[transcript]
            except KeyError:
                print('No rate4site scores are available for {}'.format(
                    transcript))
                transcript_r4s_scores = [np.nan] * len(transcript_pep)

            # calculate expected counts for each codon
            transcript_cds = transcript_cds[:-3]  # remove the stop codon
            codon_mutation_rates = seq_utils.get_codon_mutation_rates(transcript_cds)
            all_cds_ns_counts = seq_utils.count_poss_ns_variants(transcript_cds)
            all_cds_ns_sites = seq_utils.count_ns_sites(transcript_cds)

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
                total_exp_mis_counts = enst_mp_counts[transcript][-2]
            except KeyError:
                print(
                    'Transcript {} not found in {}'.format(
                        transcript, configs['enst_mp_counts'])
                )
                continue

            # permutation test
            codon_mis_probs = [x[1] for x in codon_mutation_rates]
            p = codon_mis_probs / np.sum(codon_mis_probs)
            pmt_matrix = seq_utils.permute_missense(
                total_exp_mis_counts, len(transcript_pep), p
            )

            # get the PDB ID and PDB chain associated with this transcript
            try:
                pdb_id, pdb_chain = enst_to_pdb[transcript]
            except KeyError:
                logging.critical(
                    '%s not found in given database %s',
                    transcript, configs['enst_to_pdb']
                )
                logging.critical('%s was skipped ...', transcript)
                print('%s not found in given database %s' %
                      (transcript, configs['enst_to_pdb']))
                print('%s was skipped ...' % transcript)
                continue

            # print message
            print(
                'Estimating COSMIS scores for:', transcript, ensp_id,
                 uniprot_id, pdb_id, pdb_chain
            )

            chain = pdb_utils.get_pdb_chain(
                pdb_id, pdb_chain, configs['pdb_dir']
            )

            if chain is None:
                print('No chain was available for %s' % transcript)
                continue
            
            # get all contact pairs in the PDB structure
            all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa)]
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
            pdb_to_uniprot_mapping = sifts_residue_mapping.pdb_to_uniprot(
                pdb_id, pdb_chain, uniprot_id
            )

            if uniprot_to_pdb_mapping is None or pdb_to_uniprot_mapping is None:
                print('Mapping between UniProt and PDB failed! Skipped!')
                continue

            # row identifier for each COSMIS score record
            id_fields = [transcript, ensp_id, uniprot_id]
            
            for i, a in enumerate(transcript_pep, start=1):
                # @TODO need to get the corresponding position in the PDB file
                # It is highly likely that the amino acid at position i in the
                # PDB file will not be the same amino acid at position i in the
                # ENSP peptide sequence
                try:
                    pdb_pos = uniprot_to_pdb_mapping[i]
                    res = chain[pdb_pos]
                except KeyError:
                    # mtr3d.append((i, compute_mtr1d(i, missense_counts, synonymous_counts,
                    #                            expected_counts, window=31)))
                    # print('Residue', i, 'not found in', pdb_file)
                    logging.critical(
                        'Residue %s in %s not found in chain %s in PDB file: %s', 
                        i, ensp_id, pdb_chain, pdb_id
                    )
                    continue
                    
                # check that the amino acid in ENSP sequence matches 
                # that in the PDB structure
                pdb_aa = seq1(res.get_resname())
                if a != pdb_aa:
                    logging.critical(
                        'Residue in %s and %s of %s did not match that in PDB %s '
                        'chain %s at %s: %s vs %s',
                        uniprot_id, ensp_id, transcript, pdb_id, pdb_chain, i, 
                        a, pdb_aa
                    )
                    logging.critical(
                        'Please first check UniProt sequence is identical to ENSP '
                        'sequence. If this is true, check if there is any oddity '
                        'in the SIFTS residue-level mapping.'
                    )
                    continue
                    # sys.exit(1)

                contact_res = indexed_contacts[res]
                contacts_pdb_pos = [r.get_id()[1] for r in contact_res]
                seq_seps = ';'.join(str(x) for x in
                                    [i - pdb_pos for i in contacts_pdb_pos])
                num_contacts = len(contacts_pdb_pos)

                total_missense_obs = missense_counts.setdefault(i, 0)
                total_synonymous_obs = synonymous_counts.setdefault(i, 0)
                mis_var_sites = site_variability_missense.setdefault(i, 0)
                syn_var_sites = site_variability_synonymous.setdefault(i, 0)
                try:
                    total_missense_poss = all_cds_ns_counts[i - 1][0]
                    total_synonyms_poss = all_cds_ns_counts[i - 1][1]
                    total_missense_sites = all_cds_ns_sites[i - 1][0]
                    total_synonymous_sites = all_cds_ns_sites[i - 1][1]
                    total_synonymous_rate = codon_mutation_rates[i - 1][0]
                    total_missense_rate = codon_mutation_rates[i - 1][1]
                    phylop_scores = transcript_phylop_scores[i - 1][3]
                    gerp_scores = transcript_gerp_scores[i - 1][3]
                    if not all(gerp_scores):
                        gerp_scores = [np.nan] * 3
                    r4s_score = transcript_r4s_scores[i - 1]
                except IndexError:
                    print('list index out of range:', i)
                    continue

                # get sequence context of the contact set
                try:
                    seq_context = seq_utils.get_codon_seq_context(
                        contacts_pdb_pos + [i], transcript_cds
                    )
                except IndexError:
                    break

                # compute the GC content of the sequence context
                if len(seq_context) == 0:
                    print('No nucleotides were found in sequence context!')
                    continue
                gc_fraction = seq_utils.gc_content(seq_context)

                if contacts_pdb_pos:
                    all_ensp_pos = []
                    for j in contacts_pdb_pos:
                        # @TODO need to get the corresponding position in ENSP
                        # amino acid sequence using the SIFTS mapping table
                        try:
                            ensp_pos = pdb_to_uniprot_mapping[j]
                        except KeyError:
                            logging.critical(
                                'PDB residue %s in %s chain %s not found '
                                'in %s, skipped',
                                j, pdb_id, pdb_chain, ensp_id
                            )
                            continue

                        # count the total # observed variants in contacting residues
                        total_missense_obs += missense_counts.setdefault(ensp_pos, 0)
                        total_synonymous_obs += synonymous_counts.setdefault(ensp_pos, 0)
                        mis_var_sites += site_variability_missense.setdefault(ensp_pos, 0)
                        syn_var_sites += site_variability_synonymous.setdefault(ensp_pos, 0)
                        # count the total # expected variants
                        try:
                            total_missense_poss += all_cds_ns_counts[ensp_pos - 1][0]
                            total_synonyms_poss += all_cds_ns_counts[ensp_pos - 1][1]
                            total_missense_sites += all_cds_ns_sites[ensp_pos - 1][0]
                            total_synonymous_sites += all_cds_ns_sites[ensp_pos - 1][1]
                            total_synonymous_rate += codon_mutation_rates[ensp_pos - 1][0]
                            total_missense_rate += codon_mutation_rates[ensp_pos - 1][1]
                        except IndexError:
                            logging.critical(
                                'PDB residue %s in %s chain %s out of the range '
                                'of %s codons, something must be wrong with SIFTS mapping',
                                j, pdb_id, pdb_chain, transcript
                            )
                            continue
                        all_ensp_pos.append(ensp_pos)

                    # get permutation statistics
                    pmt_mean, pmt_sd, p_value = seq_utils.get_permutation_stats(
                        pmt_matrix, all_ensp_pos + [i], total_missense_obs
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
                            '{:.3f}'.format(pmt_mean),
                            '{:.3f}'.format(pmt_sd),
                            '{:.3e}'.format(p_value),
                            '{:.3f}'.format(np.mean(phylop_scores)),
                            '{:.3f}'.format(np.mean(gerp_scores)),
                            r4s_score,
                            enst_mp_counts[transcript][2],
                            enst_mp_counts[transcript][4],
                            total_exp_mis_counts,
                            len(transcript_pep),
                        ]
                    )
                else:
                    continue

            with open(file=cosmis_file, mode='wt') as opf:
                header = get_dataset_headers()
                csv_writer = csv.writer(opf, delimiter='\t')
                csv_writer.writerow(header)
                csv_writer.writerows(cosmis)


if __name__ == '__main__':
    main()
