import filecmp
import pytest

from pathlib import Path
from hmm_from_fasta.curate_protein_sequences import *


def test_get_input_paths():
    """
    Test that if input files are found, they are all added to the list of paths
    """
    test_data_path = Path(__file__).parent.parent / "data" / "test_data"
    expected = [test_data_path / "less_than_sequence_min_length", test_data_path / "repeated_records", test_data_path / "large_data_file",
                test_data_path / "contains_phages"]

    result = get_input_paths()

    for x in range(len(expected)):
        assert expected[x] == result[x]


def test_filter_out_phages():
    """
    Test that all records that contain keyword "phage" are no longer found in data after filter_out_phages step.
    """
    result = filter_out_phages()

    for record in result:
        assert record.id != "YP_009218802.1"


def test_remove_duplicates():
    """
    Test that all duplicates are removed
    """
    record_list = filter_out_phages()
    remove_duplicates(record_list)
    result_path = Path(
        __file__).parent.parent / "data" / "results" / "temp_files" / "filtered_fasta_protein"

    record_list_repeats = {}

    with open(result_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in record_list_repeats:
                record_list_repeats[record.id] = 1
            else:
                record_list_repeats[record.id] += 1

    for repeat_value in record_list_repeats.values():
        assert repeat_value == 1


def test_min_length():
    record_list = filter_out_phages()
    remove_duplicates(record_list)
    result_path = Path(
        __file__).parent.parent / "data" / "results" / "temp_files" / "filtered_fasta_protein"

    with open(result_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            assert len(record.seq) > 1


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
