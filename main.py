from typing import Iterator, List

from Bio.SeqRecord import SeqRecord

from benchmarks.protein_match_benchmark import benchmark_protein_match
from gene_prediction.gene_predictor import get_annotated_genes
from transcription.transcriptor import transcribe_genes_into_rna
from translation.protein_translator import (
    save_protein_sequences_to_fasta,
    translate_multiple_rna_sequences,
)


def start(annotated_genes_fasta: str, whole_genome_fasta: str, protein_fasta: str):

    annotated_genes: Iterator[SeqRecord] = get_annotated_genes(
        whole_genome_fasta, annotated_genes_fasta
    )

    rna_sequences: List[SeqRecord] = transcribe_genes_into_rna(annotated_genes)

    protein_sequences: List[str] = translate_multiple_rna_sequences(rna_sequences)

    save_protein_sequences_to_fasta(
        protein_sequences, "./data/generated_data/proteins.fasta"
    )

    benchmark_protein_match(protein_sequences, protein_fasta)
