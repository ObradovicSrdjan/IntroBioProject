import argparse
from Bio import SeqIO
from cds_finder import get_coding_sequences
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
        "--protein_fasta",
        nargs="?",
        default=DEFAULT_ECOLI_PROTEIN_FASTA_FILE,
        type=str,
        help=f"Path to the protein FASTA file (default: {DEFAULT_ECOLI_PROTEIN_FASTA_FILE})",
    )

    args = parser.parse_args()

    # Process genome coding sequences
    sequence = next(SeqIO.parse(args.genome_fasta, "fasta"))
    print(f"Finding coding sequences in {args.genome_fasta}...")
    cds_list = get_coding_sequences(sequence)
    print(f"Found {len(cds_list)} coding sequences.")
    for idx, cds in enumerate(cds_list[:5]):  # Preview first 5 sequences
        print(f"CDS {idx + 1}: {cds[:30]}...")

    # Translate CDS into proteins
    proteins = translate_multiple_dna_sequences(cds_list)
    print(f"Translated {len(proteins)} proteins.")
    save_protein_sequences_to_fasta(
        proteins, file_path="generated_data/ecoli/predicted_proteins.fasta"
    )
    for idx, protein in enumerate(proteins[:5]):  # Show first 5 proteins
        print(f"Protein {idx + 1}: {protein[:30]}...")

    # Process protein sequences
    # print(f"\nParsing protein sequences from {args.protein_fasta}...")
    # protein_sequences = parse_protein_fasta(args.protein_fasta)
    # print(f"Found {len(protein_sequences)} protein sequences.")
    # for idx, protein in enumerate(protein_sequences[:5]):  # Preview first 5 sequences
    #     print(f"Protein {idx + 1}: {protein[:30]}...")

    # # Compare translated proteins with those from the FASTA file
    # matches = compare_proteins(proteins, protein_sequences)
    # print_comparison_results(matches)


if __name__ == "__main__":
    main()
