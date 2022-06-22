from cosmis.mutation_rates.trinucleotide_context_rates import MUTATION_RATES
from cosmis.utils.genetic_code import GENETIC_CODE


class ContactSet:
    """

    """
    def __init__(self, central_res, contact_res):
        """

        Parameters
        ----------
        central_res : Residue
            The central amino acid residue of the contact set.
        contact_res : list
            A list of residues that are in contact with the central residue.
        """
        self.central_res = central_res
        self.contact_res = contact_res

    def get_size(self):
        """
        Get the size of the contact set (number of residues).

        Returns
        -------
        int
            Number of residues in the contact set.
        """
        return len(self.contact_res) + 1

    def get_codons(self, enst_id):
        """

        Parameters
        ----------
        enst_id : str
            Ensembl transcript for the nucleotide sequence encoding this
            protein.

        Returns
        -------
        list
            A list containing the codons for each of the amino acid in the
            contact set.

        """
        # get the coding nucleotide sequence
        pass

    def __init__(self, cds, codon_ids, transcript_id):
        """

        Parameters
        ----------
        cds : str
            Coding sequence.
        codon_ids :
            A set of codon ID numbers.
        transcript_id : str
            Ensembl transcript identifier.
        """
        self._cds = cds
        self._codon_ids = codon_ids
        self._transcript_id = transcript_id
        self._mutation_rates = None
        self._gnomad_variants = None

    def get_mutation_rates(self):
        """

        Returns
        -------

        """
        if self._mutation_rates is None:
            synonymous_rate = 0
            nonsynonymous_rate = 0
            for i in self._codon_ids:
                x, y = self._compute_codon_mutation_rates(i)
                synonymous_rate += x
                nonsynonymous_rate += y
            self._mutation_rates = synonymous_rate, nonsynonymous_rate
        return self._mutation_rates

    def get_gnomad_variants(self, gnomad_variant_dict):
        """

        Parameters
        ----------
        gnomad_variants

        Returns
        -------

        """
        if self._gnomad_variants is None:
            # extract gnomad variants
            try:
                all_variants = gnomad_variant_dict[self._transcript_id]['variants']
            except KeyError:
                print('No variants found in %s in gnomAD', self._transcript_id)
        nonsynonymous_variants = []
        synonymous_variants = []
        for variant in all_variants:
            w = variant[0]  # wild-type amino acid
            v = variant[-1]  # mutant amino acid
            pos = int(variant[1:-1])  # position in the protein sequence
            if pos in self._codon_ids:
                if w != v:  # missense variant
                    nonsynonymous_variants.append(variant)
                else:  # synonymous variant
                    synonymous_variants.append(variant)
        return synonymous_variants, nonsynonymous_variants

    def _compute_codon_mutation_rates(self, codon_number):
        """

        Parameters
        ----------
        codon_number

        Returns
        -------

        """
        # one nucleotide before and one nucleotide after the codon
        sequence_context = self._cds[(codon_number - 1) * 3 - 1:(codon_number - 1) * 3 + 4]

        # mutate nucleotide in the codon iteratively
        synonymous_rate = 0
        nonsynonymous_rate = 0
        for i in range(3):
            trinucleotide = sequence_context[i:i + 3]
            rates = MUTATION_RATES[trinucleotide]
            for k, v in rates.items():
                if GENETIC_CODE[trinucleotide] == GENETIC_CODE[k]:
                    synonymous_rate += v
                else:
                    nonsynonymous_rate += v
        return synonymous_rate, nonsynonymous_rate

