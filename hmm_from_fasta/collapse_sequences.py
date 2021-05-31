import subprocess


from pathlib import Path


SEQUENCE_ID_THRESHOLD = 0.9
SEQUENCE_ID = 1.0
NUM_CORES = 8


def cluster_sequences(input_fasta_file) -> Path:
    """
    Takes in a fasta file, minimum fraction coverage, minimum fraction identity,calls cd-hit

    cd-hit collapses the input sequences into non-redundant representatives at the specifies levels
    
    :param input_fasta_file: filtered fasta file, no clusters
    :return: clustered_fasta_file: collapsed fasta file, clustered by cd-hit
    """
    clustered_fasta_file = Path(input_fasta_file).parent / "clustered_fasta_file"

    cd_hit_cmd = "cd-hit -i %s -o %s -c %s -s %s" % (
        input_fasta_file, clustered_fasta_file, SEQUENCE_ID_THRESHOLD, SEQUENCE_ID)
    subprocess.run(cd_hit_cmd.split())

    return clustered_fasta_file


def all_by_all_blast(clustered_fasta_file) -> Path:
    """
    Takes a clustered fasta file as input, and formats file to be a BLAST protein database

    runs BLAST on the file to itself as the BLAST-formatted database

    :param clustered_fasta_file:
    :return: tab-delimited formatted BLAST results file
    """

    format_database_cmd = "makeblastdb -in %s -dbtype prot" % clustered_fasta_file
    subprocess.run(format_database_cmd.split())

    blast_results_file = Path(clustered_fasta_file).parent / "blast_results"

    all_by_all_blast_cmd = "blastp -query %s -out %s -db %s -outfmt 6 -num_threads %s" % (
        clustered_fasta_file, blast_results_file, clustered_fasta_file, NUM_CORES)
    subprocess.run(all_by_all_blast_cmd.split())

    return blast_results_file
