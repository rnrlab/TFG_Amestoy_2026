#### Renaming Fasta and ExplorePipolin Files ####

"""
This script uses ExplorePipolin to annotate presence and structure of pipolins found in
fasta nucleotide sequences (in our case from bacterial genomes).

ExplorePipolin is meant to work with fasta files (.fna) whose file name does not contain
more than 16 characters. We need to make sure this is true for our files, specially since
NCBI-datasets output files tend to have >16 characters.

In our case, parallelization is needed since we have >16000 bacterial genomes
"""

import os
import sys
import shutil
import subprocess

#=========================================================================
# Step 0: Rename fasta files (.fna)
#=========================================================================

# Check if the correct arguments have been passed
if len(sys.argv) != 4:
    print("Usage: python ExplorePipolin.py <base_directory> <tmp_directory> <final_directory>")
    sys.exit(1)

# Main directory as second argument
base_directory = sys.argv[1]
tmp_directory = sys.argv[2]
final_directory = sys.argv[3]

# Check if the directory exists
if not os.path.exists(base_directory) or not os.path.isdir(base_directory):
    print(f"Error: The directory {base_directory} does not exist.")
    sys.exit(1)

if not os.path.exists(tmp_directory):
    os.makedirs(tmp_directory)

# Traverse directories and files within the main directory (data/)
for root, dirs, files in os.walk(base_directory):
    for file in files:
        if file.endswith(".fna"):
            file_path = os.path.join(root, file) # Coge el separador correcto
            dir_name = os.path.basename(root) # Coge el último componente de la ruta
            new_file_name = f"{dir_name}.fna"
            new_file_path = os.path.join(tmp_directory, new_file_name)
            try:
                shutil.copy2(file_path, new_file_path)
                print(f"{file_path} -> {new_file_path}")
            except Exception as e:
                print(f"{file_path} -> {e}")

#=========================================================================
# Step 1: Run ExplorePipolin
#=========================================================================

# Include '-j <int>' and '--threads <int>' to define number of processes in parallel and threads in use
EP_cmd = f'find {tmp_directory}/*.fna | parallel -j 4 --verbose "explore_pipolin {{}} --cpus 8 --keep-tmp --out-dir {final_directory}/{{/.}}"'
subprocess.run(f"mkdir -p {final_directory}", shell=True)
subprocess.run(EP_cmd, shell=True)