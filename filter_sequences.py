from Bio import SeqIO
from pathlib import Path


def filter_sequences_on_coverage(fasta_file, coverage_heuristics_dict):
    """
    Takes a fasta file and a dictionary containing coverage heuristics information

    any sequences that don't subscribe to these heuristics are removed

    :param fasta_file: fasta file containing unfiltered protein sequence records
    :param coverage_heuristics_dict: dictionary has the form {k1:v1,k2:v2,...,kn:vn}
    - kx represents 100 residues
    - vx represents the min fraction relative length between each sequence and the median sequence length of the cluster
    :return: filtered_fasta_file
    """
    record_sequence_lengths = {}
    sequence_lengths = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        record_sequence_lengths[record.id] = len(record.seq)
        sequence_lengths.append(len(record.seq))

    sequence_lengths.sort()
    median = sequence_lengths[int(len(sequence_lengths)/2)]
    coverage_key = int(median/100)

    if coverage_key in coverage_heuristics_dict:
        coverage_threshold = (1.0 + coverage_heuristics_dict[coverage_key])/2
    else:
        coverage_threshold = (1.0 + max(coverage_heuristics_dict.values()))/2

    for sequence_length in record_sequence_lengths:
        if sequence_length.value < coverage_threshold * median or sequence_length.value * coverage_threshold > median:
            record_sequence_lengths.pop(sequence_length.key)

    if len(record_sequence_lengths) == 0:
        return None
    elif len(record_sequence_lengths) == len(sequence_lengths):
        return fasta_file
    else:
        output_path = Path(fasta_file).parent / "filtered_fasta.faa"
        filtered_fasta_file = open(output_path, "w")

        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in record_sequence_lengths:
                filtered_fasta_file.write(str(record) + "\n")
        filtered_fasta_file.close()

        return filtered_fasta_file

