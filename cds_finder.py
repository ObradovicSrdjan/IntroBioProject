from Bio import SeqIO
from typing import List
from Bio.SeqRecord import SeqRecord

def get_coding_sequences(sequence: SeqRecord) -> List[str]:
    """Finds coding sequences (CDS) in a DNA sequence.

    Args:
        sequence (SeqRecord): A Biopython SeqRecord containing the DNA sequence.

    Returns:
        List[str]: List of CDS sequences from start codon (ATG) to stop codon.
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    coding_sequences = []
    start_index = 0

    while (start_index := sequence.seq.find(start_codon, start_index)) != -1:
        stop_index = -1

        for stop_codon in stop_codons:
            temp_stop_index = sequence.seq.find(stop_codon, start_index)
            if temp_stop_index != -1 and (stop_index == -1 or temp_stop_index < stop_index):
                stop_index = temp_stop_index

        if stop_index != -1 and (stop_index - start_index) > 5:
            coding_sequence = str(sequence.seq[start_index:stop_index + 3])
            if coding_sequence not in coding_sequences:
                coding_sequences.append(coding_sequence)

        start_index += len(start_codon)

    return coding_sequences
