from Bio.SeqRecord import SeqRecord


def transcribe_dna_into_rna(dna_sequence: SeqRecord) -> SeqRecord:
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
    rna_seq_record = SeqRecord(rna_sequence, id=dna_sequence.id, description="RNA sequence")

    return rna_seq_record
