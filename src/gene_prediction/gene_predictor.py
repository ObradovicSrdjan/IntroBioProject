import logging
import os
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_annotated_genes(
    whole_genome_fasta: str, annotated_genes_fasta: str
) -> Iterator[SeqRecord]:
    if os.path.exists(annotated_genes_fasta):
        logging.info(f"Loading annotated genes from: {annotated_genes_fasta}...")
        annotated_genes = SeqIO.parse(annotated_genes_fasta, "fasta")

        logging.debug(f"The first 5 gene sequences are:")
        for gene in list(annotated_genes)[:5]:
            logging.debug(f"Gene sequence: {str(gene.seq)[:20]}...")

        # Reset the iterator
        annotated_genes = SeqIO.parse(annotated_genes_fasta, "fasta")
        return annotated_genes
    else:
        logging.error(
            f"Gene prediction method not implemented yet. Please provide a valid FASTA file"
        )
        exit(1)
