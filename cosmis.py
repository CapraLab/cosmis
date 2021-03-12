#!/usr/bin/env python3

"""
To add a brief summary of this script.
"""

import os, csv
import gzip
import json
import logging
import urllib
import numpy as np
from collections import defaultdict
from argparse import ArgumentParser
from cosmis.mapping.sifts import SIFTS
from cosmis.utils import pdb_utils, seq_utils
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import is_aa

from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore', BiopythonWarning)


def compute_mtr1d(pos, ns_counts, syn_counts, expected_counts, window=31):
    """

    Parameters
    ----------
    pos : int
        Sequence position.
    ns_counts : dict

    syn_counts : dict

    expected_counts : list

    window : int
        Window size.

    Returns
    -------
    float


    """
    total_ns_obs = 0
    total_syn_obs = 0
    total_ns_exp = 0
    total_syn_exp = 0
    # handle the special case of the first 15 residues
    if pos < 16:
        for i in range(31):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            total_ns_exp += expected_counts[i][0]
            total_syn_exp += expected_counts[i][1]
    # handle the special case of the last 15 residues
    elif pos > len(expected_counts) - 15:
        for i in range(len(expected_counts) - 31, len(expected_counts)):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            total_ns_exp += expected_counts[i][0]
            total_syn_exp += expected_counts[i][1]
    else:
        for i in range(pos - window // 2, pos + window // 2 + 1):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            if i > 1 or i < len(expected_counts):
                total_ns_exp += expected_counts[i - 1][0]
                total_syn_exp += expected_counts[i - 1][1]
    try:
        mtr1d = (total_ns_obs / (total_ns_obs + total_syn_obs)) / \
            (total_ns_exp / (total_ns_exp + total_syn_exp))
    except ZeroDivisionError:
        mtr1d = 0

    return total_ns_obs, total_syn_obs, mtr1d


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
        '--mtr1d', dest='mtr1d', required=False, default=False, action='store_true', 
        help='''If specified, computes COSMIS scores (a.k.a. MTR) implemented 
        according to Traynelis et al., Genome Research, 2017.'''
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
        vv, ac, an = variant
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        pos = vv[1:-1]  # position in the protein sequence
        # only consider rare variants
        if int(ac) / int(an) > 0.001:
            continue
        if w != v:  # missense variant
            missense_counts[int(pos)] += 1
        else:  # synonymous variant
            synonymous_counts[int(pos)] += 1
    return missense_counts, synonymous_counts


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

    # get the directory where all output files will be stored
    output_dir = os.path.abspath(configs['output_dir'])

    # create SIFTS mapping table
    sifts_residue_mapping = SIFTS(configs['sifts_residue_mapping_file'])

    # set output file suffix
    if args.mtr1d:
        suffix = '_mtr1ds.tsv'
    else:
        suffix = '_cosmis.tsv'

    # compute the COSMIS scores for each transcript
    with open(args.transcripts, 'rt') as ipf:
        for transcript in ipf:
            transcript = transcript.strip()
            print('Processing transcript %s' % transcript)
            cosmis = []
            cosmis_file = os.path.join(output_dir, transcript + suffix)
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
            if len(transcript_cds) / 3 != len(transcript_pep) + 1:
                logging.critical(
                    'Incomplete CDS for', transcript, '. Skipped.'
                )
                continue

            # check that the CDS does not contain invalid nucleotides
            if any([x not in {'A', 'T', 'C', 'G'} for x in set(transcript_cds)]):
                logging.critical('Invalid CDS! Skipped.')
                print(transcript_cds)
                continue

            # calculate expected counts for each codon
            codon_mutation_rates = seq_utils.get_codon_mutation_rates(transcript_cds)
            all_cds_ns_counts = seq_utils.count_cds_ns(transcript_cds)

            # tabulate variants at each site
            missense_counts, synonymous_counts = count_variants(variants)

            # compute the total number of missense variants
            total_mis_counts = 0
            for k, v in missense_counts.items():
                total_mis_counts += v

            # permutation test
            permuted_missense_mutations = seq_utils.permute_missense(
                total_mis_counts, len(transcript_pep)
            )

            # only compute MTR1D scores if asked on the command-line
            if args.mtr1d:
                """
                @TODO to set the number of expected counts here
                """
                expected_counts = None
                for i, a in enumerate(transcript_pep, start=1):
                    cosmis.append(
                        [transcript, ensp_id, uniprot_id, i, a] +
                        list(
                            compute_mtr1d(
                                i, missense_counts, synonymous_counts,
                                expected_counts, window=31
                            )
                        )
                    )
                    
                with open(cosmis_file, mode='wt') as opf:
                    header = [
                        'TRANSCRIPT', 'PROTEIN', 'UNIPROT', 'ENSP_POS', 
                        'ENSP_AA', 'MISSENSE', 'SYNONYMOUS', 'MTR1D'
                    ]
                    csv_writer = csv.writer(opf, delimiter='\t')
                    csv_writer.writerow(header)
                    csv_writer.writerows(cosmis)
                print('Finished calcualting MTR1D scores for %s!' % transcript)
                continue

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
                    cosmis.append(id_fields + [i, a] + [np.nan] * 10)
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
                    cosmis.append(
                        id_fields + [i, a, pdb_pos, pdb_aa] + [np.nan] * 8
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
                try:
                    total_missense_poss = all_cds_ns_counts[i - 1][0]
                    total_synonyms_poss = all_cds_ns_counts[i - 1][1]
                    total_synonymous_rate = codon_mutation_rates[i - 1][0]
                    total_missense_rate = codon_mutation_rates[i - 1][1]
                except IndexError:
                    print('list index out of range:', i)
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
                        # count the total # expected variants
                        try:
                            total_missense_poss += all_cds_ns_counts[ensp_pos - 1][0]
                            total_synonyms_poss += all_cds_ns_counts[ensp_pos - 1][1]
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
                    contact_res_indices = [pos - 1 for pos in all_ensp_pos] + [i - 1]
                    mean_missense_counts = np.mean(
                        permuted_missense_mutations[:, contact_res_indices].sum(axis=1)
                    )
                    sd_missense_counts = np.std(
                        permuted_missense_mutations[:, contact_res_indices].sum(axis=1)
                    )

                    # push results for the current residue
                    cosmis.append(
                        id_fields + 
                        [i, a, pdb_pos, pdb_aa] + 
                        [pdb_id, pdb_chain] + 
                        [
                            seq_seps, 
                            num_contacts,
                            total_synonyms_poss,
                            total_missense_poss,
                            '{:.3e}'.format(total_synonymous_rate),
                            total_synonymous_obs, 
                            '{:.3e}'.format(total_missense_rate),
                            total_missense_obs,
                            '{:.3f}'.format(mean_missense_counts),
                            '{:.3f}'.format(sd_missense_counts)
                        ]
                    )
                else:
                    continue

            with open(file=cosmis_file, mode='wt') as opf:
                header = [
                    'enst_id', 'ensp_id', 'uniprot_id', 'ensp_pos', 'ensp_aa', 
                    'pdb_pos', 'pdb_aa',  'pdb_id', 'chain_id', 'seq_separations',
                    'num_contacts', 'synonymous_poss', 'missense_poss',
                    'synonymous_rate', 'synonymous_obs', 'missense_rate',
                    'missense_obs', 'permutation_mean', 'permutation_sd'
                ]
                csv_writer = csv.writer(opf, delimiter='\t')
                csv_writer.writerow(header)
                csv_writer.writerows(cosmis)


if __name__ == '__main__':
    main()
