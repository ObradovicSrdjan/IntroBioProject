import os
from typing import List
from Bio.SeqRecord import SeqRecord


codon_to_amino_acid = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",  # Isoleucine
    "ATG": "M",  # Methionine (Start codon)
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",  # Threonine
    "AAC": "N",
    "AAT": "N",  # Asparagine
    "AAA": "K",
    "AAG": "K",  # Lysine
    "AGC": "S",
    "AGT": "S",  # Serine
    "AGA": "R",
    "AGG": "R",  # Arginine
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",  # Leucine
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",  # Proline
    "CAC": "H",
    "CAT": "H",  # Histidine
    "CAA": "Q",
    "CAG": "Q",  # Glutamine
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",  # Arginine
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",  # Valine
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",  # Alanine
    "GAC": "D",
    "GAT": "D",  # Aspartic Acid
    "GAA": "E",
    "GAG": "E",  # Glutamic Acid
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",  # Glycine
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",  # Serine
    "TTC": "F",
    "TTT": "F",  # Phenylalanine
    "TTA": "L",
    "TTG": "L",  # Leucine
    "TAC": "Y",
    "TAT": "Y",  # Tyrosine
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",  # Stop codons
    "TGC": "C",
    "TGT": "C",  # Cysteine
    "TGG": "W",  # Tryptophan
}


def translate_rna_to_protein(rna_sequence: str):
    protein_sequence = ""
    for i in range(0, len(rna_sequence) - 2, 1):
        codon = rna_sequence[i : i + 3]
        if codon in codon_to_amino_acid:
            if codon_to_amino_acid[codon] == "M":
                protein_sequence, end_of_transcription = read_orfs(rna_sequence[i:])
    return protein_sequence


def read_orfs(sliced_dna_sequence: str):
    protein_sequence = ""
    for i in range(0, len(sliced_dna_sequence) - 2, 3):
        codon = sliced_dna_sequence[i : i + 3]
        if codon in codon_to_amino_acid:
            if codon_to_amino_acid[codon] == "*":  # Stop translation at stop codon
                return protein_sequence, i + 3
            protein_sequence += codon_to_amino_acid[codon]
        else:
            protein_sequence += "X"  # Unknown codon

    return protein_sequence + "!", -1


def translate_multiple_rna_sequences(coding_sequences: List[SeqRecord]):
    return [
        translate_rna_to_protein(str(sequence.seq)) for sequence in coding_sequences
    ]


def save_protein_sequences_to_fasta(protein_sequences, file_path="proteins.fasta"):
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w") as fasta_file:
        for i, protein in enumerate(protein_sequences):
            fasta_file.write(f">Protein_{i+1}\n")
            fasta_file.write(f"{protein}\n")
