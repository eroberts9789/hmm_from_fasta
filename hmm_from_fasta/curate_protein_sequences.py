import os
import sys


from Bio import SeqIO
from pathlib import Path


INPUT_PATH = Path(__file__).parent.parent / "data" / "test_data"
OUTPUT_FILE_NAME = "filtered_fasta_protein"
OUTPUT_PATH = INPUT_PATH.parent / "results" / "temp_files" / OUTPUT_FILE_NAME
SEQUENCE_MIN_LENGTH = 1


def get_input_paths() -> list:
    """
    Takes in a path to directory containing input fasta files, returns a list of paths to the fasta files

    :return: input_paths, list of paths to input protein files if any files are found
    """
    input_files = os.listdir(INPUT_PATH)
    input_paths = [INPUT_PATH / file for file in input_files]

    if input_paths:
        return input_paths
    else:
        print("Error: no files found in input directory")
        sys.exit(1)


def filter_out_phages() -> list:
    """
    Parses through records in input files, appends records to filtered_sequences if keyword "phage" not found in record

    :return: filtered_records, list of filtered records
    """
    input_paths = get_input_paths()
    filtered_records = list()

    for input_path in input_paths:
        with open(input_path) as handle:

            for record in SeqIO.parse(handle, "fasta"):
                if "phage" not in record.description:
                    filtered_records.append(record)
    return filtered_records


def remove_duplicates(filtered_records) -> Path:
    """
    Removes duplicates in filtered_records list, writes all records in list to OUTPUT_PATH

    verifies length of each sequence is longer than SEQUENCE_MIN_LENGTH

    :param filtered_records: list of records from all protein files without keyword "phage"
    :return: Path to filtered file
    """
    record_ids = []
    records_to_output = []

    for record in filtered_records:
        if record.id not in record_ids and len(record.seq) > SEQUENCE_MIN_LENGTH:
            records_to_output.append(record)
            record_ids.append(record.id)

    SeqIO.write(records_to_output, OUTPUT_PATH, "fasta")

    return OUTPUT_PATH







