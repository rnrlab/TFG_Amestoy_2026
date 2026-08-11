#### Identification of integrons within bacterial genomes containing piPolBs ####

"""
We start with the .fna nucleotide FASTA files located within the ncbi_dataset/data/GCA* folders.
In our case, the path is pipol_ncbi_dataset/ncbi_dataset/data/GCA*.
The files are named after the GCA assembly of their genome.

We use IntegronFinder2.0 as tool for integrons identification in bacterial genomes.
"""

import os
import subprocess

#=========================================================================
# Step 0: Stalish paths and input files
#=========================================================================

base_path = "./pipol_ncbi_dataset/ncbi_dataset/data" # Directory with genomes
base_out_dir = "./results_integronfinder"

# Create output directory if necessary
os.makedirs(base_out_dir, exist_ok=True)

#=========================================================================
# Step 1: Execute IntegronFinder2.0
#=========================================================================

integronfinder_cmd = f'''
find {base_path} -name "*.fna" | \
parallel -j 8 --verbose '
integron_finder {{}} \
--local-max \
--cpu 4 \
--evalue-attc 5 \
--outdir {base_out_dir}
'
'''
subprocess.run(integronfinder_cmd, shell=True)
