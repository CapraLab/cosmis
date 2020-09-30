#!/usr/bin/env python3

from cosmis.utils.genetic_code import GENETIC_CODE
from cosmis.mutation_rates.trinucleotide_context_rates import  MUTATION_RATES_UNIQUE


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
            elif new_aa != aa and new_aa != 'Stop':
                missense += 1
    return missense, synonymous


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
            raise IndexError('Codon number %s out of range.' % codon_nums)
    else:
        for i in codon_nums:
            if i > codon_nums:
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
            codon_nums += cds[((i - 1) * 3 - 1):(i * 3) + 1]
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
