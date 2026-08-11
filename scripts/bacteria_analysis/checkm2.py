#### Completeness and Contamination Level Annotation with CheckM2 ####

"""
This script uses nucleotide sequences in FASTA format (.fna files) as input
for CheckM2.
This tool annotates completeness and contamination level of these genomes, grouped
in batches. For each batch CheckM2 gives a TSV. Then we combine all of these into
one unique TSV file.

For parallelization, we use the Pool class from the multiprocessing module.
You could also try GNU parallel, but it did not work for me :]
"""

import sys
import os
import subprocess
import shutil
from multiprocessing import Pool
import glob
import pandas as pd

#=========================================================================
# Step 0: Creating variables
#=========================================================================
# Error checking for command-line use
if len(sys.argv) != 3:
    print("Usage: python checkm2.py <directory_with_genomes> <output_directory>")
    sys.exit(1)

# Variables
base_directory = sys.argv[1]
output_base = sys.argv[2]
batch_size = 10
threads = 8
parallel_batches = 4

#=========================================================================
# Step 1: Creating batches of genomes
#=========================================================================

# Listing of subdirectories with genomes
all_genomes = [os.path.join(base_directory, d) for d in os.listdir(base_directory) if d.startswith("GCA")]
print(f"-> Genomes found: {len(all_genomes)}")

# Creating batches
batches = []
for idx in range(0, len(all_genomes), batch_size):
    batches.append(all_genomes[idx:idx+batch_size])
print(f"-> Batches: {len(batches)}")
print(f"-> Batch-size: {batch_size}")

#=========================================================================
# Step 2: CheckM2 processing batches function
#=========================================================================

def process_batch(args):
    """
    Process a batch of genomes by running CheckM2 on them.

    This function receives a batch index and a list of genome directories.
    For each batch, it:
    - Creates an output directory for CheckM2 results.
    - Creates a temporary directory containing symbolic links to `.fna` files
      found in each genome directory.
    - Executes CheckM2 using the configured number of threads.
    - Removes the temporary directory once processing is complete.

    Parameters:
        args (tuple): A tuple containing:
            - i (int): Batch index.
            - batch (list of str): List of paths to genome directories, each expected to contain at least one `.fna` file.

    Returns:
        None: The function produces output files in the batch-specific output directory but does not return any value.
    """
    # args is a tuple with: i = batch number, batch = list of genomes
    i, batch = args

    # We create an output folder, format --> checkm2_out/batch_1/
    output_dir = os.path.join(output_base, f"batch_{i}")
    os.makedirs(output_dir, exist_ok=True)

    # We create the temporary folder that will serve as input for CheckM2 --> temp_batch_1/
    batch_dir = f"temp_batch_{i}"
    os.makedirs(batch_dir, exist_ok=True)

    # Creation of symlinks --> link that points to the genome directory
    for d in batch: # We traverse the genomes in the batch (list of genomes), where d is the path to a genome
        # We look for the .fna file within the directory
        files = [f for f in os.listdir(d) if f.endswith(".fna")]
        if not files:
            print(f"No .fna files in {d}")
            continue

        genome_file = os.path.join(d, files[0]) # We take the path of the .fna file
        symlink_path = os.path.join(batch_dir, os.path.basename(genome_file)) # Path to symlink
        if not os.path.exists(symlink_path): # We check that the path does not exist
            os.symlink(os.path.abspath(genome_file), symlink_path) # We don't copy the genome, we only create a link to it

    # Execute CheckM2
    command =["checkm2", "predict", "--input", batch_dir, "--output-directory", output_dir, "--threads", str(threads), "-x", "fna"]
    print(f"-> Batch {i+1}...")
    subprocess.run(command, check=True)
    print(f"-> Batch {i+1} done")

    # Delete temporary folder from batch
    if os.path.exists(batch_dir):
        shutil.rmtree(batch_dir)


#=========================================================================
# Step 3: Parallel batch processing with multiprocessing
#=========================================================================

with Pool(parallel_batches) as p:
    p.map(process_batch, enumerate(batches))
    # enumerate() adds an index to each batch --> necessary since process_batch function takes a tuple: (idx, batch)
    # p.map() executes parallel batches: process_batch((i, batch))

#=========================================================================
# Step 4: Combine results
#=========================================================================

all_results = []
for batch_direc in glob.glob(os.path.join(output_base, "batch_*")): # glob.glob() returns a list with all the directories matching the pattern
    result_file = os.path.join(batch_direc, "quality_report.tsv") # Open CheckM2 results file of each batch
    if os.path.exists(result_file):
        df = pd.read_csv(result_file, sep = "\t") # Read the results file and turn into DataFrame
        all_results.append(df) # Add DataFrame to a list of results

if all_results:
    final_df = pd.concat(all_results, ignore_index=True) # Merge DataFrames
    final_file = os.path.join(output_base, "checkm2_all_genomes.tsv")
    final_df.to_csv(final_file, sep = "\t", index=False) # Save combined dataframe
    print(f"-> Results saved to {final_file}")
else:
    print(f"-> No results found")