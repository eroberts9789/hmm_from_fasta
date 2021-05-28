import subprocess


from pathlib import Path


INFLATION_NUM = None


def blast_to_mcl(all_by_all_blast_results_file, polyprotein_sequences):
    """
    Converts sequences not included in polyprotein_sequences to a .abc file

    calls mcxload on .abc file to generate a .mci and .tab file

    calls mcl on .tab file to generate newline-separated clusters

    :param all_by_all_blast_results_file:blast file produced in all_by_all blast step
    :param polyprotein_sequences: list of polyprotein like sequences to not be included in output
    :return: mcl_file_path to file containing newline-separated clusters
    """
    abc_file_path = Path(all_by_all_blast_results_file).parent / "abc_file.abc"
    mci_file_path = Path(all_by_all_blast_results_file).parent / "mci_file.mci"
    tab_file_path = Path(all_by_all_blast_results_file).parent / "tab_file.tab"
    mcl_file_path = Path(all_by_all_blast_results_file).parent / "mcl_file.mcl"

    blast_file = open(all_by_all_blast_results_file)
    abc_file = open(abc_file_path, 'w')

    for line in blast_file:
        data = line.rstrip().split()
        query = data[0]
        subject = data[1]
        expected_value = data[10]

        if query not in polyprotein_sequences and subject not in polyprotein_sequences:
            abc_line = '\t'.join([query, subject, expected_value]) + "\n"
            abc_file.write(abc_line)

    abc_file.close()

    mcxload_cmd = "mcxload --stream-mirror -abc %s -o %s -write-tab %s" \
                  % (abc_file_path, mci_file_path, tab_file_path)
    subprocess.run(mcxload_cmd.split())

    if not INFLATION_NUM:
        inflation_args = ""
    else:
        inflation_args = "-I %s " % INFLATION_NUM

    mcl_cmd = "mcl %s -use-tab %s %s -o %s" % (mci_file_path, tab_file_path, inflation_args, mcl_file_path)
    subprocess.run(mcl_cmd.split())

    return mcl_file_path


def mcl_to_fasta(mcl_file_path):
    pass
