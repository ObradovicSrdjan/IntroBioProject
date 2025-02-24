import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from cds_finder import get_coding_sequences
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
        "--protein_fasta",
        nargs="?",
        default=DEFAULT_ECOLI_PROTEIN_FASTA_FILE,
        type=str,
        help=f"Path to the protein FASTA file (default: {DEFAULT_ECOLI_PROTEIN_FASTA_FILE})",
    )

    args = parser.parse_args()

    if os.path.exists(RNA_FILE):
        print(f"Loading RNA sequence from: {RNA_FILE}...")
        whole_genome_rna_sequence: SeqRecord = next(SeqIO.parse(RNA_FILE, "fasta"))
    else:
        whole_genome_dna_sequence: SeqRecord = next(
            SeqIO.parse(args.genome_fasta, "fasta")
        )
        print(f"Transcribing DNA into RNA for: {args.genome_fasta}...")
        whole_genome_rna_sequence: SeqRecord = transcribe_dna_into_rna(
            dna_sequence=whole_genome_dna_sequence
        )
        SeqIO.write(whole_genome_rna_sequence, RNA_FILE, "fasta")

    print(f"Finding coding sequences...")
    cds_list = get_coding_sequences(whole_genome_rna_sequence)

    # Translate CDS into proteins
    proteins = translate_multiple_dna_sequences(cds_list)
    print(f"Translated {len(proteins)} proteins.")
    save_protein_sequences_to_fasta(
        proteins, file_path="generated_data/ecoli/predicted_proteins.fasta"
    )
    for idx, protein in enumerate(proteins[:5]):  # Show first 5 proteins
        print(f"Protein {idx + 1}: {protein[:30]}...")


if __name__ == "__main__":
    main()
