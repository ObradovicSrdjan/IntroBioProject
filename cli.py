import argparse
import logging

from pipeline_runner import run


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline to do translation/transcription."
    )
    # parser.add_argument(
    #     "--genome_fasta",
    #     nargs="?",
    #     default=DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE,
    #     type=str,
    #     help=f"Path to the genome FASTA file (default: {DEFAULT_ECOLI_WHOLE_GENOME_FASTA_FILE})",
    # )

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

    run()


if __name__ == "__main__":
    main()
