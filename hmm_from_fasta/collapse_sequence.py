import subprocess


from pathlib import Path


output_path = Path(__file__).parent.parent / "data"
SEQUENCE_ID_THRESHOLD = 0.9
SEQUENCE_ID = 1.0


def collapse_sequence(input_fasta_file, collapsed_fasta_file):
    """
    Takes in a fasta file, minimum fraction coverage, minimum fraction identity,calls cd-hit

    cd-hit collapses the input sequences into non-redundant representatives at the specifies levels

    :return: collapsed fasta file grouped by clusters
    """
    collapsed_fasta_file = output_path / collapsed_fasta_file
    cd_hit_cmd = "cd-hit -i %s -o %s -c %s -G %s" % (input_fasta_file, collapsed_fasta_file, SEQUENCE_ID_THRESHOLD, SEQUENCE_ID)
    subprocess.run(cd_hit_cmd.split())
    return collapsed_fasta_file
