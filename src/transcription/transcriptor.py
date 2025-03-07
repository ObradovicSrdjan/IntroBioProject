import logging
from typing import Iterator, List

from Bio.SeqRecord import SeqRecord


def transcribe_dna_sequence_into_rna(dna_sequence: SeqRecord) -> SeqRecord:
    """
    Translates a DNA sequence into an RNA sequence.
    We assume that the input DNA sequence is the coding strand.
    The RNA sequence is transcribed from the template strand,
    thus producing a sequence with uracil (U) instead of thymine (T).

    Args:
        sequence (SeqRecord): A SeqRecord object containing a DNA sequence.

    Returns:
        SeqRecord: A SeqRecord object containing the corresponding RNA sequence.
    """
    if dna_sequence.seq is None:
        raise ValueError("The sequence attribute is None.")

    # Transcribe the DNA sequence to RNA
    rna_sequence = dna_sequence.seq.replace("T", "U")

    # Create a new SeqRecord for the RNA sequence
    rna_seq_record = SeqRecord(
        rna_sequence, id=dna_sequence.id, description="RNA sequence"
    )

    return rna_seq_record


def transcribe_genes_into_rna(genes: Iterator[SeqRecord]) -> List[SeqRecord]:
    """
    Transcribes a list of DNA sequences into RNA sequences.

    Args:
        genes (List[SeqRecord]): A list of SeqRecord objects containing DNA sequences.

    Returns:
        List[SeqRecord]: A list of SeqRecord objects containing the corresponding RNA sequences.
    """
    logging.info("Transcribing DNA into RNA...")
    rna_sequences = [transcribe_dna_sequence_into_rna(gene) for gene in genes]
    logging.debug(f"The first 5 RNA sequences are:")
    for rna in rna_sequences[:5]:
        logging.debug(f"RNA sequence: {str(rna.seq)[:20]}...")
    return rna_sequences
