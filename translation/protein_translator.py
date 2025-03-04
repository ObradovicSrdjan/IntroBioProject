import os
from typing import List
from Bio.SeqRecord import SeqRecord
import logging


codon_to_amino_acid = {
    "AUA": "I",
    "AUC": "I",
    "AUU": "I",  # Isoleucine
    "AUG": "M",  # Methionine (Start codon)
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",  # Threonine
    "AAC": "N",
    "AAU": "N",  # Asparagine
    "AAA": "K",
    "AAG": "K",  # Lysine
    "AGC": "S",
    "AGU": "S",  # Serine
    "AGA": "R",
    "AGG": "R",  # Arginine
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",  # Leucine
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",  # Proline
    "CAC": "H",
    "CAU": "H",  # Histidine
    "CAA": "Q",
    "CAG": "Q",  # Glutamine
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",  # Arginine
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",  # Valine
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",  # Alanine
    "GAC": "D",
    "GAU": "D",  # Aspartic Acid
    "GAA": "E",
    "GAG": "E",  # Glutamic Acid
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",  # Glycine
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",  # Serine
    "UUC": "F",
    "UUU": "F",  # Phenylalanine
    "UUA": "L",
    "UUG": "L",  # Leucine
    "UAC": "Y",
    "UAU": "Y",  # Tyrosine
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",  # Stop codons
    "UGC": "C",
    "UGU": "C",  # Cysteine
    "UGG": "W",  # Tryptophan
}


def translate_rna_to_protein(rna_sequence: str) -> str:
    protein_sequence = read_orfs(rna_sequence)

    # set first char to be M
    protein_sequence = "M" + protein_sequence[1:]

    # Remove stop codons
    protein_sequence = protein_sequence.replace("*", "")

    return protein_sequence


def read_orfs(rna_sequence: str):
    protein_sequence = ""
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i : i + 3]
        if codon in codon_to_amino_acid:
            protein_sequence += codon_to_amino_acid[codon]
        else:
            protein_sequence += "X"  # Unknown codon

    return protein_sequence


def translate_multiple_rna_sequences(coding_sequences: List[SeqRecord]) -> List[str]:
    logging.info("Translating RNA into proteins...")
    protein_sequences = [
        translate_rna_to_protein(str(sequence.seq)) for sequence in coding_sequences
    ]
    logging.debug(f"The first 5 protein sequences are:")
    for protein in protein_sequences[:5]:
        logging.debug(f"Protein sequence: {protein}")
    return protein_sequences


def save_protein_sequences_to_fasta(
    protein_sequences, file_path="proteins.fasta"
) -> None:
    logging.info(f"Saving protein sequences to {file_path}")
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w") as fasta_file:
        for i, protein in enumerate(protein_sequences):
            fasta_file.write(f">Protein_{i+1}\n")
            fasta_file.write(f"{protein}\n")
    logging.info(f"Proteins saved to {file_path}")
