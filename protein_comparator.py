from Bio import pairwise2


def compare_proteins(translated_proteins, fasta_proteins):
    """
    Compare translated proteins from CDS with proteins from the FASTA file.
    Returns a list of matches with alignment scores.
    """
    matches = []
    for idx, translated_protein in enumerate(translated_proteins):
        best_match = None
        best_score = 0
        for fasta_protein in fasta_proteins:
            alignments = pairwise2.align.globalxx(
                translated_protein, fasta_protein, one_alignment_only=True
            )
            if alignments:
                score = alignments[0].score
                if score > best_score:
                    best_score = score
                    best_match = fasta_protein
        if best_match:
            matches.append((translated_protein, best_match, best_score))
    return matches


def print_comparison_results(matches):
    """Prints the best matches with scores."""
    print("\nProtein Comparison Results:")
    for idx, (translated, fasta, score) in enumerate(matches[:5]):  # Show top 5 matches
        print(f"Match {idx + 1}: Score {score}")
        print(f"Translated: {translated[:30]}...")
        print(f"Fasta: {fasta[:30]}...\n")
