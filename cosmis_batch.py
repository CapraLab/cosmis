#!/usr/bin/env python3

import sys, os, csv
import gzip
import json
import urllib
from collections import defaultdict
from argparse import ArgumentParser
from mtr3d.struct.contact import Contact, BACKBONE_ATOMS
from mtr3d.mapping.sifts import SIFTS
from mtr3d.mapping.ensembl_uniprot_pdb import EnsemblUniProtPDB
from mtr3d.utils import pdb_utils
from mtr3d.utils.genetic_code import GENETIC_CODE
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, NeighborSearch, is_aa
from mtr3d.mutation_rates.trinucleotide_context_rates import MUTATION_RATES_UNIQUE


def get_codon_mutation_rates(cds):
    """

    Parameters
    ----------
    cds : str
        Coding sequence.

    Returns
    -------
    mutation probabilities : list
        A list of tuples consisting of synonymous and nonsynonymous mutation
        probabilities.
    """
    if len(cds) % 3 != 0:
        raise ValueError('Given CDS length is not a multiple of 3.')

    num_codons = len(cds) // 3

    mutation_rates = []
    # one nucleotide before and one nucleotide after the codon
    for codon_number in range(1, num_codons + 1):
        # determine the codon sequence
        codon_sequence = cds[(codon_number - 1) * 3:codon_number * 3]

        synonymous_rate = 0
        nonsynonymous_rate = 0

        # determine the mutation rate of the first and the last codons
        # consider only two mutatable nucleotides
        if codon_number == 1 or codon_number == num_codons:
            # first codon
            if codon_number == 1:
                sequence_context = cds[:4]
                # i is the zero-indexed position of the mutated nucleotide
                for i in range(1, 3):
                    trinucleotide = sequence_context[i - 1:i + 2]
                    rates = MUTATION_RATES_UNIQUE[trinucleotide]
                    for k, v in rates.items():
                        mutant_sequence = codon_sequence[:i] + k[1] + codon_sequence[i + 1:]
                        if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                            synonymous_rate += v
                        elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                            and GENETIC_CODE[mutant_sequence] != 'STOP':
                            nonsynonymous_rate += v
            # last codon
            else:
                sequence_context = cds[-4:]
                # i is the zero-indexed position of the mutated nucleotide
                for i in range(0, 2):
                    trinucleotide = sequence_context[i:i + 3]
                    rates = MUTATION_RATES_UNIQUE[trinucleotide]
                    for k, v in rates.items():
                        mutant_sequence = codon_sequence[:i] + k[1] + codon_sequence[i + 1:]
                        if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                            synonymous_rate += v
                        elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                            and GENETIC_CODE[mutant_sequence] != 'STOP':
                            nonsynonymous_rate += v

        # codons other than the first and the last
        # consider all three mutatable nucleotides
        else:
            # one nucleotide before and one nucleotide after the codon
            sequence_context = cds[(codon_number - 1) * 3 - 1:(codon_number - 1) * 3 + 4]
            # mutate nucleotide in the codon iteratively
            for i in range(3):
                trinucleotide = sequence_context[i:i + 3]
                rates = MUTATION_RATES_UNIQUE[trinucleotide]
                for k, v in rates.items():
                    codon_seq_list = list(codon_sequence)
                    codon_seq_list[i] = k[1]
                    mutant_sequence = ''.join(codon_seq_list)
                    if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                        synonymous_rate += v
                    elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                        and GENETIC_CODE[mutant_sequence] != 'STOP':
                        nonsynonymous_rate += v

        mutation_rates.append((synonymous_rate, nonsynonymous_rate))
    return mutation_rates


def search_for_all_contacts(residues, radius=6):
    """
    Search for all contacts in the given set of residues based on
    distances between heavy atoms (atoms other than hydrogen atoms).

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
    all_contacts = [Contact(res_a=c[0], res_b=c[1]) for c in ns.search_all(radius, level='R')]
    return all_contacts


def parse_cmd():
    """

    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-i', '--input', dest='input', required=True,
                        type=str, help='''Ensembl transcript ID of the 
                        transcript for which to compute a MTR3D profile.''')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='''Whether to overwrite 
                        already computed MTR3D scores.''')
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true', help='''Whether to output verbose
                        data: number of contacting residues and number of 
                        missense and synonymous variants in the neighborhood
                        of the mutation site.''')
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
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


def get_uniprot_aa_seq(uniprot_id):
    """

    Parameters
    ----------
    uniprot_id

    Returns
    -------

    """
    from urllib.request import HTTPError
    fasta_url = 'https://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
    try:
        url_stream = urllib.request.urlopen(fasta_url)
    except HTTPError:
        return None
    aa_seq = ''
    for l in url_stream.read().decode('ascii').split('\n'):
        if not l.startswith('>'):
            aa_seq += l.strip()
    return aa_seq


def get_pdb_chain(pdb_file, pdb_chain):
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
    # read in the PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    try:
        structure = pdb_parser.get_structure(id="given_pdb", file=pdb_file)
    except (FileNotFoundError, ValueError) as e:
        print('PDB file cannot be retrieved', pdb_file)
        return None
    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        print('No chain ' + pdb_chain + ' was found in ' + pdb_file)
        return None
    return chain


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
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd()

    # parse configuration file
    configs = parse_config(args.config)

    # ENSEMBL cds
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

    # directory where to store the output files
    output_dir = os.path.abspath(configs['output_dir'])

    # read transcript to swiss model mapping
    with open(args.input, 'rt') as ipf:
        transcripts = [l.strip().split() for l in ipf]
        
    # compute the MRT scores
    for x, y in transcripts:
        transcript = x
        pdb_file = 'data/swiss-model/SWISS-MODEL_Repository/' + y + '.pdb'
        pdb_chain = y[-1]
    
        if os.path.exists(
            os.path.join(output_dir, transcript + '_mtr3ds.tsv')
        ) and not args.overwrite:
            continue

        print('Processing transcript %s' % transcript)

        mtr3ds = []

        # get the amino acid sequence of the transcript
        try:
            # Ensembl peptide ID for the transcript
            ensp_id = transcript_variants[transcript]['ensp'][0]
        except KeyError:
            print(
                'Transcript %s not found in %s' % 
                (transcript, configs['gnomad_variants'])
            )
            print('Skip to the next transcript...')
            continue

        # get the peptide sequence from peptide sequence database
        transcript_pep_seq = get_transcript_pep_seq(
            transcript, ensp_id, pep_dict
        )

        # get all variants of this transcript reported in gnomAD
        try:
            variants = transcript_variants[transcript]['variants']
        except KeyError:
            print('No variants found for %s in gnomAD', transcript)
            print('Skip to the next transcript...')
            continue

        # get the coding sequence of the transcript
        try:
            transcript_cds = ensembl_cds_dict[transcript].seq
        except KeyError:
            print(
                '''No CDS found in Ensembl CDS database! 
                Looking for it in the CCDS database ...'''
            )
            transcript_cds = None

        if transcript_cds is None:
            try:
                ccds_id = transcript_variants[transcript]['ccds'][0]
                transcript_cds = ccds_dict[ccds_id].seq
            except KeyError:
                print('ERROR: No CDS found in CCDS database!')
                print('Skip to the next transcript...')
                continue

        if transcript_pep_seq is None:
            print('ERROR: No peptide sequence found in Ensembl database!')
            print('Skip to the next transcript...')
            continue

        # stop if the CDS is invalid
        if len(transcript_cds) / 3 != len(transcript_pep_seq) + 1:
            print('ERROR: incomplete CDS for', transcript)
            print('Skip to the next transcript...')
            continue

        # check that the CDS does not contain invalid nucleotides
        if any([x not in {'A', 'T', 'C', 'G'} for x in set(transcript_cds)]):
            print('ERROR: invalid CDS for', transcript_cds)
            print('Skip to the next transcript...')
            continue

        # print message
        print(
            'Estimating MRT3D scores for:', 
            transcript, 
            ensp_id,
            pdb_file
        )

        chain = get_pdb_chain(pdb_file, pdb_chain)
        
        if chain is None:
            print('ERROR: %s not found in structure: %s!' % (pdb_chain, pdb_file))
            print('Skip to the next transcript...')
            continue

        all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa)]
        all_contacts = search_for_all_contacts(all_aa_residues, radius=8)
        
        # calculate expected counts for each codon
        codon_mutation_rates = get_codon_mutation_rates(transcript_cds)

        if len(codon_mutation_rates) < len(all_aa_residues):
            print('ERROR: peptide sequence has less residues than structure!')
            print('Skip to the next transcript...')
            continue
        
        # tabulate variants at each site
        # missense_counts and synonymous_counts are dictionary that maps
        # amino acid positions to variant counts
        missense_counts, synonymous_counts = count_variants(variants)

        # index all contacts by residue ID
        indexed_contacts = defaultdict(list)
        for c in all_contacts:
            indexed_contacts[c.get_res_a()].append(c.get_res_b())
            indexed_contacts[c.get_res_b()].append(c.get_res_a())

        valid_case = True
        for seq_pos, seq_aa in enumerate(transcript_pep_seq, start=1):
            # check that the amino acid in ENSP sequence matches 
            # that in the PDB structure
            try:
                res = chain[seq_pos]
            except KeyError:
                print('PDB file is missing residue:', seq_aa, 'at', seq_pos)
                continue
            pdb_aa = seq1(res.get_resname())
            if seq_aa != pdb_aa:
                print('Residue in ENSP did not match that in PDB at', seq_pos)
                print('Skip to the next transcript...')
                valid_case = False
                break

            contact_res = indexed_contacts[res]
            contacts_pdb_pos = [r.get_id()[1] for r in contact_res]
            
            seq_seps = ';'.join(
                str(x) for x in [i - seq_pos for i in contacts_pdb_pos]
            )

            total_missense_obs = missense_counts.setdefault(seq_pos, 0)
            total_synonymous_obs = synonymous_counts.setdefault(seq_pos, 0)
            total_synonymous_rate = codon_mutation_rates[seq_pos - 1][0]
            total_missense_rate = codon_mutation_rates[seq_pos - 1][1]
            for j in contacts_pdb_pos:
                # count the total # observed variants in contacting residues
                total_missense_obs += missense_counts.setdefault(j, 0)
                total_synonymous_obs += synonymous_counts.setdefault(j, 0)

                # count the total # expected variants
                try:
                    total_synonymous_rate += codon_mutation_rates[j - 1][0]
                    total_missense_rate += codon_mutation_rates[j - 1][1]
                except IndexError:
                    valid_case = False
                    break

            # compute the fraction of expected missense variants
            mtr3ds.append(
                [
                    transcript, ensp_id, seq_pos, seq_aa, seq_seps, 
                    '{:.2e}'.format(total_synonymous_rate), 
                    total_synonymous_obs,
                    '{:.2e}'.format(total_missense_rate), 
                    total_missense_obs
                ]
            )

        if not valid_case:
            continue

        with open(
                file=os.path.join(output_dir, transcript + '_mtr3ds.tsv'), 
                mode='wt'
            ) as opf:
            header = [
                'transcript_id', 'peptide_id', 'position', 'amino_acid', 
                'contacts', 'synonymous_rate', 'synonymous_count', 
                'missense_rate', 'missense_count'
            ]
            csv_writer = csv.writer(opf, delimiter='\t')
            csv_writer.writerow(header)
            csv_writer.writerows(mtr3ds)


if __name__ == '__main__':
    main()
