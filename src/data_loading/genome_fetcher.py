from typing import List
from Bio import Entrez
import logging
from Bio.Entrez.Parser import DictionaryElement
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = (
    "Your.Name.Here@example.org"  # Why do we need to provide an email address?
)


MAX_GENOMES = 2

FILE_FORMAT = "fasta"

RAW_DATA_FOLDER_PATH = "data/raw"


def get_genome_ids(organism: str, max_results: int = MAX_GENOMES) -> List[str]:
    search_query = f'"{organism}"[Organism] AND "complete genome"[Title]'
    handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    if isinstance(record, DictionaryElement) and "IdList" in record:
        id_list = record["IdList"]
        logging.info(f"Found {len(id_list)} genomes for {organism}")
        if id_list:
            return record["IdList"]
    logging.error(f"No genome found for {organism}")
    return []


def download_genomes(organisms: List[str], file_format: str = FILE_FORMAT) -> None:
    logging.info("Dowloading genes from Entrez...")
    for organism in organisms:
        genome_ids = get_genome_ids(organism)
        for genome_id in genome_ids:
            handle = Entrez.efetch(
                db="nucleotide", id=genome_id, rettype=FILE_FORMAT, retmode="text"
            )
            filename = f"{RAW_DATA_FOLDER_PATH}/genome_{genome_id}.{FILE_FORMAT}"
            with open(filename, "w") as f:
                f.write(handle.read())
            logging.info(f"Genome saved to {filename}")
            handle.close()


def load_genomes(organisms: List[str]) -> List[SeqRecord]:
    # Todo fix this
    # if not os.path.exists(RAW_DATA_FOLDER_PATH) or not any(
    #     os.scandir(RAW_DATA_FOLDER_PATH)
    # ):
    #     download_genomes(organisms)
    # else:
    #     logging.info("Genomes already downloaded, loading from disk...")

    download_genomes(organisms)

    genomes = []
    for file in os.listdir(RAW_DATA_FOLDER_PATH):
        if file.endswith(FILE_FORMAT):
            genomes.append(SeqIO.read(f"{RAW_DATA_FOLDER_PATH}/{file}", FILE_FORMAT))

    return genomes
