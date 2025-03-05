from typing import Iterator, List, Tuple

from Bio.SeqRecord import SeqRecord

START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}


def get_coding_sequences(genes: Iterator[SeqRecord]) -> List[str]:
    coding_sequences = []
    for gene in genes:
        coding_sequences += get_coding_sequences_for_gene(gene)
    return coding_sequences


def get_coding_sequences_for_gene(gene: SeqRecord) -> List[str]:
    coding_sequences = []
    whole_genome_sequence: str = str(gene.seq)

    coding_sequence = ""
    end_of_cd = -1
    for i in range(0, len(whole_genome_sequence) - 2, 1):
        codon = whole_genome_sequence[i : i + 3]
        if codon == START_CODON:
            coding_sequence, end_of_cd = read_orf(whole_genome_sequence[i:])
            coding_sequences.append(coding_sequence)

        i = end_of_cd  # Fast forward the lookup
    return coding_sequences


def read_orf(sliced_dna_sequence: str) -> Tuple[str, int]:

    if sliced_dna_sequence[0:3] != START_CODON:
        raise ValueError("No start codon found")

    coding_sequence = ""

    for i in range(0, len(sliced_dna_sequence) - 2, 3):
        codon = sliced_dna_sequence[i : i + 3]
        coding_sequence += codon
        if codon in STOP_CODONS:
            return coding_sequence, i + 3

    return coding_sequence + "!", -1
