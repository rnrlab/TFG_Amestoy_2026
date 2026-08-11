#### Identifying Conjugative Elements with MacSyFinder-CONJScan ####

"""
We start with the protein fasta files (FAA) within the ncbi_dataset/data/GCA* folders.
These files have all the same name, protein.faa, but they are in different
subdirectories named after their respective genome assembly.

This FAA files serve as input for the macsyfinder model, CONJScan.
"""

import subprocess

#=========================================================================
# Step 0: Stablish paths
#=========================================================================

base_path = "pipol_ncbi_dataset/ncbi_dataset/data"
conjscan_model_path = "macsyfinder_models"
base_out_dir = "results_conjscan"

#=========================================================================
# Step 1: Execute CONJScan
#=========================================================================

# Include 'parallel -j <int>' and 'macsyfinder --threads <int>' to define number of processes in parallel and threads in use
conjscan_cmd = f'''
find {base_path} -name "protein.faa" | \
parallel --verbose '
genome=$(basename $(dirname {{}}))
macsyfinder \
--sequence-db {{}} \
--models CONJScan \
--models-dir {conjscan_model_path} \
--out-dir {base_out_dir}/$genome \
--db-type ordered_replicon
'
'''
subprocess.run(conjscan_cmd, shell=True)