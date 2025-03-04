import logging
from Bio import SeqIO
from typing import Iterator, List
from Bio.SeqRecord import SeqRecord
from transcription.transcriptor import transcribe_genes_into_rna
from translation.protein_translator import (
    translate_multiple_rna_sequences,
    save_protein_sequences_to_fasta,
)
from gene_prediction.gene_predictor import get_annotated_genes


def start(annotated_genes_fasta: str, whole_genome_fasta: str, protein_fasta: str):
    """
    First step is finding the segments of the whole genome that get transcribed into RNA. (???)
    We will then find the coding sequences in the RNA. (exons/introns)
    Then we can translate the coding sequences into proteins. (amino acid table)
    Then we can compare the proteins with the proteins from the protein FASTA file.
    """

    annotated_genes: Iterator[SeqRecord] = get_annotated_genes(annotated_genes_fasta)

    rna_sequences: List[SeqRecord] = transcribe_genes_into_rna(annotated_genes)

    protein_sequences: List[str] = translate_multiple_rna_sequences(rna_sequences)

    save_protein_sequences_to_fasta(
        protein_sequences, "./generated_data/proteins.fasta"
    )

    # Compare the proteins with the proteins from the protein FASTA file
    logging.info(
        "Comparing the proteins with the proteins from the protein FASTA file..."
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
