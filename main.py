import os
import logging
from Bio import SeqIO
from typing import Iterator, List
from Bio.SeqRecord import SeqRecord
from coding_sequence_finder import get_coding_sequences
from transcription.transcriptor import transcribe_dna_into_rna
from protein_parser import parse_protein_fasta
from protein_translator import (
    translate_multiple_rna_sequences,
    save_protein_sequences_to_fasta,
)


def start(annotated_genes_fasta: str, whole_genome_fasta: str):
    """
    First step is finding the segments of the whole genome that get transcribed into RNA. (???)
    We will then find the coding sequences in the RNA. (exons/introns)
    Then we can translate the coding sequences into proteins. (amino acid table)
    Then we can compare the proteins with the proteins from the protein FASTA file.
    """

    if os.path.exists(annotated_genes_fasta):
        logging.info(f"Loading annotated genes from: {annotated_genes_fasta}...")
        annotated_genes: Iterator[SeqRecord] = SeqIO.parse(
            annotated_genes_fasta, "fasta"
        )
    else:
        logging.error(f"Gene prediction method not implemented yet.")
        exit(1)

    for gene in annotated_genes:
        logging.info(f"Gene: {str(gene.seq)[:20]}...")

    logging.info("Transcribing DNA into RNA...")
    rna_sequences: List[SeqRecord] = [
        transcribe_dna_into_rna(gene) for gene in annotated_genes
    ]
    for rna in rna_sequences[:5]:
        # Log the first 5 RNA sequences only the first 20 characters
        logging.info(f"RNA sequence: {str(rna.seq)[:20]}...")

    logging.info("Translating RNA into proteins...")
    protein_sequences = translate_multiple_rna_sequences(rna_sequences)
    for protein in protein_sequences[:5]:
        logging.info(f"Protein sequence: {protein}")

    # Save the protein sequences to a FASTA file
    # save_protein_sequences_to_fasta(protein_sequences)
