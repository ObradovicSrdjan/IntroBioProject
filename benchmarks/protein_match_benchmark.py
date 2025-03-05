import logging
from typing import List

from Bio import SeqIO


def benchmark_protein_match(protein_sequences: List[str], protein_fasta: str):

    logging.info(
        "Comparing the predicted proteins with the proteins from the protein FASTA file..."
    )
    annotated_proteins = SeqIO.parse(protein_fasta, "fasta")

    matches = 0
    non_matches = 0
    for protein in annotated_proteins:
        if str(protein.seq) in protein_sequences:
            matches += 1
            logging.info(f"Protein: {protein.id} MATCHED")
        else:
            non_matches += 1
            logging.info(f"Protein: {protein.id} NOT MATCHED")

    logging.info(
        f"matches = {matches} not matches = {non_matches} total = {matches + non_matches}"
    )
