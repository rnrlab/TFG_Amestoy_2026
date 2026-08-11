#### Script to run PanTA add over bacterial genomes ####

"""
PanTA can simply run using GFF3 (--gff) as input files (this must your usual GFF3 from prodigal + the nucleotide FASTA
sequence at the bottom, with the header "##FASTA".
Nevertheless, when you have a ton of genomes like I do PanTA crushes, so you must give not the GFF3 files alone,
but a TSV file (--tsv). This is a unique file, that for each row it has the following structure:
    GCA_xxxxx  GCA_xxxxx.gff   GCA_xxxxx.fna
Yes, for some reason now it requires the FNA file.

Here we put the batches manually 'cause we only have 2.

P.D. DO NOT SPLIT PARALOGS if you value your time and hard drive
"""

import subprocess

#=========================================================================
# Step 0 - Stablish paths and functions
#=========================================================================
dest_dir = "results_panta"
temp_dir = "temp_panta"

def run_cmd(cmd, log_file):
    """
    This function runs the cmd that you put in it
    Then returns the LOG
    :param cmd: command to run
    :param log_file: path to log file
    :return: results of cmd and log file
    """
    full_cmd = f"set -o pipefail; {cmd} 2>&1 | tee -a {log_file}"
    result = subprocess.run(full_cmd, shell=True, executable='/bin/bash')

    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")

#=========================================================================
# Step 1 - Run panta main
#=========================================================================
run_cmd(
    f"panta main --outdir {dest_dir} --tsv {temp_dir}/initial_batch.tsv "
    "--identity 0.70 --AL 0.70 --AS 0.70 "
    "--dont-split --threads 32",
    f"{temp_dir}/initial_batch_restart.log"
)

#=========================================================================
# Step 2 - Run panta add
#=========================================================================
run_cmd(
    f"panta add --collection-dir {dest_dir} --tsv {temp_dir}/batch_1.tsv "
    "--identity 0.70 --AL 0.70 --AS 0.70 "
    "--dont-split --threads 32",
    f"{temp_dir}/batch_1.log"
)