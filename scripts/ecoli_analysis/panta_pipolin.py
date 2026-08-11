##### Extraction of pipolin's GFF from ExplorePipolin and PanTA #####

import os
import shutil
import glob
import subprocess

#=========================================================================
# Step 0 - Input and Output paths
#=========================================================================
# Copying
ep_dir = "/mnt/disk2/pipolin_screening_jorge/ecoli_results_EP"
pipolin_genomes_file = "ecoli_pipolin_genomes.txt"

gff_dir = "pipolin_gff_for_panta"
os.makedirs(gff_dir, exist_ok=True)

# Pangenomeing
panta_dest_dir = "results_panta_pipolin"
os.makedirs(panta_dest_dir, exist_ok=True)

temp_dir = "temp_panta"
os.makedirs(temp_dir, exist_ok=True)

"""
#=========================================================================
# Step 1 - Copy GFFs
#=========================================================================
# Open file with pipolin-containing genomes
with open(pipolin_genomes_file, 'r') as f:
    pipolin_genomes = [line.strip().split()[0] for line in f if line.strip()]

# Copy
copied = 0
for genome in pipolin_genomes:
    pattern = f"{ep_dir}/{genome}/{genome}/pipolins/{genome}_0v0*.gff"
    files = glob.glob(pattern)

    if not files:
        print(f"Could not find GFF file for {genome} in {ep_dir}")
        continue

    source = files[0]
    dest = f"{gff_dir}/{genome}.gff"

    shutil.copy(source, dest)
    copied += 1

print(f"\n Copied {copied} GFFs")
"""

#=========================================================================
# Step 2 - Run PanTA
#=========================================================================
def run_cmd(cmd, log_file):
    full_cmd = f"set -o pipefail; {cmd} 2>&1 | tee -a {log_file}"
    result = subprocess.run(full_cmd, shell=True, executable='/bin/bash')

    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")

gff_files = glob.glob(f"{gff_dir}/*.gff")
if not gff_files:
    raise RuntimeError(f"Could not find a single GFF file in {gff_dir}")

# Construct list of GFFs
gff_list = " ".join(gff_files)

# Run PanTA main
run_cmd(
    f"panta main --outdir {panta_dest_dir} --gff {gff_list} "
    "--identity 0.70 --AL 0.70 --AS 0.70 "
    "--dont-split --threads 32 "
    "--core 0.99 --soft 0.95 --shell 0.01 ",
    f"{temp_dir}/pipolin_panta.log"
)

