from Bio import SeqIO
from typing import List

def parse_protein_fasta(fasta_path: str) -> List[str]:
    """
    Parses a protein FASTA file and returns a list of protein sequences as strings.

    :param fasta_path: Path to the protein FASTA file.
    :return: List of protein sequences.
    """
    protein_sequences = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]
    return protein_sequences