#!/usr/bin/env python3

from cosmis.utils.genetic_code import GENETIC_CODE
from cosmis.mutation_rates.trinucleotide_context_rates import MUTATION_RATES_UNIQUE
import numpy as np


def place_holder():
    """

    Returns
    -------

    """
    pass


def main():
    """

    Returns
    -------

    """
    pass


def count_codon_ns(codon):
    """
    Count the respective number of possible missense and synonymous SNPs of the
    given codon.

    Parameters
    ----------
    codon : str
        A three-letter genetic codon.

    Returns
    -------
    tuple
        A pair of # missense variants vs # synonymous variants.

    """
    nucleotides = {'A', 'T', 'C', 'G'}
    aa = GENETIC_CODE[codon]
    missense = 0
    synonymous = 0
    for i, x in enumerate(codon):
        for n in nucleotides - {x}:
            new_codon = codon[:i] + n + codon[i+1:]
            new_aa = GENETIC_CODE[new_codon]
            if new_aa == aa:
                synonymous += 1
            elif new_aa != aa and new_aa != '*':
                missense += 1
    return missense, synonymous


def is_valid_cds(seq_record):
    """
    Determines whether a given sequence is a valid transcript sequence or not.
    A valid transcript sequence is defined as one that starts with ATG, ends
    with a stop codon, that whose length is a multiple of three.

    Parameters
    ----------
    seq_record : SeqRecord
        An object of type SeqRecord.

    Returns
    -------
    bool
        True if the given sequence is a valid transcript sequence else False.

    """
    # stop codons
    stop_codons = {'TAG', 'TAA', 'TGA'}

    # the following conditionals skip noncanonical transcripts
    id_fields = seq_record.id.split('|')
    if not seq_record.seq[:3] == 'ATG':
        print('The CDS of', id_fields[0], 'does not start with ATG.')
        return False
    if seq_record.seq[-3:] not in stop_codons:
        print('The CDS of', id_fields[0], 'does not have a stop codon.')
        return False
    if len(seq_record.seq) % 3 != 0:
        print('The CDS of', id_fields[0], 'is not a multiple of 3.')
        return False

    return True


def count_cds_ns(cds):
    """

    Parameters
    ----------
    cds : str
        Coding sequence.

    Returns
    -------
    list
        A pair of # missense variants vs # synonymous variants for each
        codon in the coding sequence.

    """
    cds_ns = []
    for i in range(0, len(cds), 3):
        cds_ns.append(count_codon_ns(cds[i:i+3]))
    return cds_ns


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


def get_codon_seq_context(codon_nums, cds):
    """

    Parameters
    ----------
    codon_nums : list/int
        A list of ints, each represents a codon number.
    cds : str
        Coding sequence

    Returns
    -------
    str
        The sequence context of codons concatenated into a string.

    """
    # make sure the codon_nums are valid
    total_num_codons = len(cds) // 3
    if isinstance(codon_nums, int):
        if codon_nums > total_num_codons:
            print('The given coding sequence has %s codons:' % total_num_codons)
            print(cds)
            raise IndexError('Codon number %s out of range.' % codon_nums)
    else:
        for i in codon_nums:
            if i > total_num_codons:
                print('The given coding sequence has %s codons:' % total_num_codons)
                print(cds)
                raise IndexError('Codon number %s out of range.' % i)

    # if only given a single codon number
    if isinstance(codon_nums, int):
        if codon_nums == 1:
            # wrap to the last nucleotide
            return cds[-1:] + cds[:4]
        elif codon_nums == total_num_codons:
            # wrap to the first nucleotide
            return cds[-4:] + cds[:1]
        else:
            return cds[((codon_nums - 1) * 3 - 1):(codon_nums * 3) + 1]
    else:
        codon_seq_context = ''
        for i in codon_nums:
            codon_seq_context += cds[((i - 1) * 3 - 1):(i * 3) + 1]
        return codon_seq_context


def gc_content(seq):
    """
    Compute the GC content of the given sequence.

    Parameters
    ----------
    seq : str
        A valid nucleotide (DNA) sequence.

    Returns
    -------
    float
        The fraction of nucleotides that are G and C.

    """
    # make sure that the input is valid
    if set(seq).difference(set('ATCG')):
        raise ValueError('Invalid input sequence: %s.' % seq)

    gc_count = 0
    for x in seq:
        if x in set('GC'):
            gc_count += 1
    return gc_count / len(seq)


def count_cg_gc(seq):
    """
    Count the number of GC/CG 2-mers in the given sequence.

    Parameters
    ----------
    seq : str
        A valid nucleotide (DNA) sequence.

    Returns
    -------
    tuple
        The number of GC/CG 2-mers in the given sequence.

    """
    cg_count = 0
    gc_count = 0
    for i in range(len(seq) - 1):
        two_mer = seq[i:(i + 2)]
        if two_mer == 'CG':
            cg_count += 1
        elif two_mer == 'GC':
            gc_count += 1
    return cg_count, gc_count


def translate(cds):
    """
    Translate the given coding sequence to amino acid sequence, dropping the
    stop codon.

    Parameters
    ----------
    cds : str
        Coding sequence

    Returns
    -------
    str
        Amino acid sequence.

    """
    if len(cds) % 3 != 0:
        raise ValueError('Invalid CDS: CDS length is not a multiple of 3.')

    aa_seq = []
    # skip the top codon
    for x in range(0, len(cds) - 3, 3):
        codon = cds[x:(x+3)].seq
        aa_seq.append(GENETIC_CODE[codon])
    return ''.join(aa_seq)


def snp_dms(cds):
    """
    Deep mutational scan of single nucleotide substitutions at each site.

    Parameters
    ----------
    cds : str
        Coding sequence.

    Returns
    -------
    list
        A list of lists where each nested list contains all the possible
        missense variants due to single nucleotide substitutions at the
        corresponding site.
    """
    # make sure the input CDS is valid
    if set(cds).difference(set('ATCG')):
        raise ValueError(
            'Invalid DNA sequence: containing letters other than ATCG.'
        )
    if len(cds) % 3 != 0:
        raise ValueError('Invalid CDS: CDS length is not a multiple of 3.')

    # now do a deep mutational scan
    dms_variants = []
    for i in range(0, len(cds) - 3, 3):
        current_codon = cds[i:i+3]
        wt_aa = GENETIC_CODE[current_codon]
        current_vars = set()
        for k, x in enumerate(current_codon):
            for y in {'A', 'T', 'C', 'G'} - {x}:
                var_codon = current_codon[:k] + y + current_codon[k+1:]
                var_aa = GENETIC_CODE[var_codon]
                # missense variant
                if var_aa != wt_aa and var_aa != '*':
                    current_vars.add(var_aa)
        dms_variants.append(sorted(current_vars))

    return dms_variants


def count_ns_sites(cds):
    """
    Count the number of nonsynonymous and synonymous sites in the given cds.

    Parameters
    ----------
    cds : str
        Coding DNA sequence.

    Returns
    -------
    tuple
        Number of nonsynonymous and synonymous sites as a tuple.

    """
    # make sure that the given CDS is valid
    if set(cds).difference(set('ATCG')):
        raise ValueError(
            'Invalid DNA sequence: containing letters other than ATCG.'
        )
    if len(cds) % 3 != 0:
        raise ValueError('Invalid CDS: CDS length is not a multiple of 3.')

    # now compute dN and dS
    n = 0
    s = 0
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        wt_aa = GENETIC_CODE[codon]
        for k, x in enumerate(codon):
            for y in {'A', 'T', 'C', 'G'} - {x}:
                var_codon = codon[:k] + y + codon[k+1:]
                var_aa = GENETIC_CODE[var_codon]
                if wt_aa != var_aa:
                    n += 1.0 / 3
                else:
                    s += 1.0 / 3
    return n, s


def get_uniprot_aa_seq(uniprot_id):
    """
    Retrieve amino acid sequence given its UniProt identifier.

    Parameters
    ----------
    uniprot_id : str
        UniProt identifier of the amino acid sequence.

    Returns
    -------
    str
        The amino acid sequence associated with the UniProt identifier.

    """
    import urllib
    from urllib.error import HTTPError
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


def get_missense_sites(variants):
    """

    Parameters
    ----------
    variants

    Returns
    -------

    """
    mis_sites = set()
    for variant in variants:
        vv, _, _ = variant
        w = vv[0]
        v = vv[-1]
        if w != v:
            mis_sites.add(int(vv[1:-1]))
    return mis_sites


def permute_missense(m, l, p=None, n=10000):
    """
    To be added ...

    Parameters
    ----------
    m
    l
    p
    n

    Returns
    -------

    """
    if not (p is None) and l != len(p):
        raise ValueError('')
    # normalize p
    if p is not None:
        p = p / np.sum(p)
    missense_matrix = []
    for _ in range(n):
        v = np.zeros(l)
        m_sites = np.random.choice(range(l), m, replace=True, p=p)
        v[m_sites] = 1
        missense_matrix.append(v)
    return np.stack(missense_matrix)

