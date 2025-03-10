import time
from typing import List
from Bio import Entrez
import os
from Bio.SeqRecord import SeqRecord
import logging
from Bio import SeqIO
from Bio.Entrez.Parser import DictionaryElement

from src.config import ENTREZ_EMAIL, FILE_FORMAT, PROTEINS_RAW_DATA_FILE_PATH


Entrez.email = ENTREZ_EMAIL


RETRY_COUNT = 5
RETRY_DELAY_SECONDS = 3


def get_protein_ids_for_accession(accession: str) -> List[str]:
    handle = Entrez.esearch(db="protein", term=accession)
    record = Entrez.read(handle)
    handle.close()
    if isinstance(record, DictionaryElement) and "IdList" in record:
        id_list = record["IdList"]
        logging.info(f"Found {len(id_list)} proteins for accession {accession}")
        if len(id_list) == 0:
            return []
        if id_list:
            return record["IdList"]
    logging.error(f"No proteins found for accession {accession}")
    return []


def download_proteins(accession) -> None:
    logging.info("Downloading proteins from Entrez...")
    for accession in accession:
        filename = f"{PROTEINS_RAW_DATA_FILE_PATH}/proteins_{accession}.{FILE_FORMAT}"
        if os.path.exists(filename):
            logging.info(
                f"Proteins for accession {accession} already exist, skipping download."
            )
            continue

        get_protein_ids = get_protein_ids_for_accession(accession)

        for attempt in range(RETRY_COUNT):
            try:
                handle = Entrez.efetch(
                    db="protein",
                    id=",".join(get_protein_ids),
                    rettype=FILE_FORMAT,
                    retmode="text",
                    retmax=2,
                )
                with open(filename, "w") as f:
                    f.write(handle.read())
                logging.info(f"Proteins saved to {filename}")
                handle.close()
                break
            except Exception as e:
                logging.error(
                    f"Failed to download proteins for accession {accession} on attempt {attempt + 1}: {e}"
                )
                if attempt < RETRY_COUNT - 1:
                    time.sleep(RETRY_DELAY_SECONDS)
                else:
                    logging.error(
                        f"Exceeded maximum retries for proteins for accession {accession}"
                    )


def get_genome_accessions(whole_genomes: List[SeqRecord]) -> List[str]:
    return [genome.description.split()[0] for genome in whole_genomes]


def load_proteins_for_genomes(
    whole_genomes: List[SeqRecord], load_from_disk: bool = False
) -> dict:
    proteins = {}

    if not load_from_disk:
        accessions: List[str] = get_genome_accessions(whole_genomes)
        download_proteins(accessions)

    for file in os.listdir(PROTEINS_RAW_DATA_FILE_PATH):
        if file.endswith(FILE_FORMAT):
            accession = file.split("_")[1].split(".")[0]
            logging.info(f"Loading proteins for {accession} from {file}")
            # Do we care about the accession for each protein?
            proteins[accession] = list(
                SeqIO.parse(f"{PROTEINS_RAW_DATA_FILE_PATH}/{file}", FILE_FORMAT)
            )

    return proteins
