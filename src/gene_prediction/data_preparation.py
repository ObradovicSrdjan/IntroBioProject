import os
import random
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

CSV_FILE_PATH = "data/generated_data/genes_training_set.csv"
NUMBER_OF_GENES = 4639
GENOME_SIZE_BP = 4641652
LONGEST_GENE_BP = 8622
MAX_START_INDEX = GENOME_SIZE_BP - LONGEST_GENE_BP - 1


def write_sequences_to_csv(is_gene: bool, sequences: List[str]):
    if not os.path.exists(CSV_FILE_PATH):
        with open(CSV_FILE_PATH, "w") as f:
            f.write("IsGene,Sequence\n")
            for sequence in sequences:
                f.write(f"{int(is_gene)},{sequence}\n")
    else:
        with open(CSV_FILE_PATH, "a") as f:
            for sequence in sequences:
                f.write(f"{int(is_gene)},{sequence}\n")


def write_annotated_genes_to_csv() -> List[str]:
    annotated_genes = SeqIO.parse(
        "data/verified_data/ecoli/genes/MG1655_annotated_genes.fna", "fasta"
    )
    gene_sequences = [gene.seq for gene in annotated_genes]
    write_sequences_to_csv(True, gene_sequences)
    return gene_sequences


def generate_non_genes(gene_sequences: List[str]) -> List[str]:
    random_non_genes = []

    for _ in range(NUMBER_OF_GENES):
        random_start_index = random.randint(0, MAX_START_INDEX)

        random_gene_length_index = random.randint(0, NUMBER_OF_GENES - 1)

        random_length = len(gene_sequences[random_gene_length_index])

        whole_genome = next(
            SeqIO.parse(
                "data/verified_data/ecoli/genome/GCF_000005845.2_ASM584v2_genomic.fna",
                "fasta",
            )
        ).seq

        random_non_gene = whole_genome[
            random_start_index : random_start_index + random_length
        ]
        if random_non_gene not in gene_sequences:
            random_non_genes.append(random_non_gene)

    print(len(random_non_genes))
    return random_non_genes


def write_non_genes_to_csv(gene_sequences: List[str]):

    random_non_genes = generate_non_genes(gene_sequences)
    write_sequences_to_csv(False, random_non_genes)


def get_details_from_csv():
    df = pd.read_csv(CSV_FILE_PATH)

    # Calculate the length of each sequence
    df["Sequence_Length"] = df["Sequence"].apply(len)

    # Plot the distribution of sequence lengths
    plt.figure(figsize=(10, 6))
    plt.hist(df["Sequence_Length"], bins=50, edgecolor="black")
    plt.title("Distribution of Sequence Lengths")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.show()

    # Plot the count of gene vs. non-gene sequences
    plt.figure(figsize=(10, 6))
    df["IsGene"].value_counts().plot(kind="bar", color=["blue", "orange"])
    plt.title("Count of Gene vs. Non-Gene Sequences")
    plt.xlabel("IsGene")
    plt.ylabel("Count")
    plt.xticks(ticks=[0, 1], labels=["Non-Gene", "Gene"], rotation=0)
    plt.show()

    # Display basic statistics about the sequences
    print(df.describe())


gene_sequences = write_annotated_genes_to_csv()
write_non_genes_to_csv(gene_sequences)
get_details_from_csv()
