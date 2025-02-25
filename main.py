import argparse
import os
import logging
from Bio import SeqIO
from typing import Iterator, List
from Bio.SeqRecord import SeqRecord
from coding_sequence_finder import get_coding_sequences
from transcriptor import transcribe_dna_into_rna
from protein_parser import parse_protein_fasta
from protein_translator import (
    translate_multiple_dna_sequences,
    save_protein_sequences_to_fasta,
)
from protein_comparator import compare_proteins, print_comparison_results


DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE = (
    "data/ecoli/genome/GCF_000005845.2_ASM584v2_genomic.fna"
)
DEFAULT_ECOLI_PROTEIN_FASTA_FILE = "data/ecoli/proteins/uniprotkb_Escherichia_coli_str_K_12_sub_2025_02_12.fasta"  # cspell:disable-line
DEFAULT_ECOLI_ANNOTATED_GENES_FASTA_FILE = "data/ecoli/genes/MG1655_annotated_genes.fna"
RNA_FILE = "generated_data/ecoli/whole_genome_rna.fasta"


def main():
    parser = argparse.ArgumentParser(
        description="Find coding and protein sequences in an E. coli genome/protein FASTA file."
    )
    parser.add_argument(
        "--genome_fasta",
        nargs="?",
        default=DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE,
        type=str,
        help=f"Path to the genome FASTA file (default: {DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE})",
    )
    parser.add_argument(
        "--annotated_genes_fasta",
        nargs="?",
        default=DEFAULT_ECOLI_ANNOTATED_GENES_FASTA_FILE,
        type=str,
        help=f"Path to the annotated genes FASTA file (default: {DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE})",
    )
    parser.add_argument(
        "--protein_fasta",
        nargs="?",
        default=DEFAULT_ECOLI_PROTEIN_FASTA_FILE,
        type=str,
        help=f"Path to the protein FASTA file (default: {DEFAULT_ECOLI_PROTEIN_FASTA_FILE})",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Increase output verbosity",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    """
        First step is finding the segments of the whole genome that get transcribed into RNA. (???)
        We will then find the coding sequences in the RNA. (exons/introns)
        Then we can translate the coding sequences into proteins. (amino acid table)
        Then we can compare the proteins with the proteins from the protein FASTA file.
    """

    if os.path.exists(args.annotated_genes_fasta):
        logging.info(f"Loading annotated genes from: {args.annotated_genes_fasta}...")
        annotated_genes: Iterator[SeqRecord] = SeqIO.parse(RNA_FILE, "fasta")
    else:
        # Do ML to predict genes
        pass

    logging.info(f"Finding coding sequences...")
    coding_sequences: List[str] = get_coding_sequences(annotated_genes)

    logging.info(f"Found {len(coding_sequences)} coding sequences.")
    for i, coding_sequence in enumerate(coding_sequences[:10]):
        logging.info(f"  {i + 1}: {coding_sequence}")

    # Save coding sequences to fasta

    # proteins = translate_coding_sequences(coding_sequences)
    # logging.info(f"Translated {len(proteins)} proteins.")
    # save_protein_sequences_to_fasta(
    #     proteins, file_path="generated_data/ecoli/predicted_proteins.fasta"
    # )


if __name__ == "__main__":
    main()
