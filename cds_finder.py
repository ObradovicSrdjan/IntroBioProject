from typing import List
from Bio.SeqRecord import SeqRecord


def get_coding_sequences(sequence: SeqRecord) -> List[str]:
    start_codon = "ATG"
    coding_sequences = []
    whole_genome_sequence: str = str(sequence.seq)

    coding_sequence = ""
    for i in range(0, len(whole_genome_sequence) - 2, 1):
        codon = whole_genome_sequence[i : i + 3]
        if codon == start_codon:
            coding_sequence, end_of_cd = read_orf(whole_genome_sequence[i:])
            coding_sequences.append(coding_sequence)

        i = end_of_cd  # Fast forward the lookup
    return coding_sequences


def read_orf(sliced_dna_sequence: str):

    stop_codons = {"TAA", "TAG", "TGA"}

    coding_sequence = ""

    for i in range(0, len(sliced_dna_sequence) - 2, 3):
        codon = sliced_dna_sequence[i : i + 3]
        coding_sequence += codon
        if codon in stop_codons:
            return coding_sequence, i + 3

    return coding_sequence + "!", -1
