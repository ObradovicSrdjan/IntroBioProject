import argparse
import logging
from main import start


DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE = (
    "data/ecoli/genome/GCF_000005845.2_ASM584v2_genomic.fna"
)
DEFAULT_ECOLI_PROTEIN_FASTA_FILE = "data/ecoli/proteins/uniprotkb_Escherichia_coli_str_K_12_sub_2025_02_12.fasta"  # cspell:disable-line
DEFAULT_ECOLI_ANNOTATED_GENES_FASTA_FILE = "data/ecoli/genes/MG1655_annotated_genes.fna"
RNA_FILE = "generated_data/ecoli/whole_genome_rna.fasta"


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline to do translation/transcription."
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

    start(args.annotated_genes_fasta, args.genome_fasta)


if __name__ == "__main__":
    main()
