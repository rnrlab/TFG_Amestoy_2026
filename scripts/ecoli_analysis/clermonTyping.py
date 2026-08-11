#### Identification of Escherichia coli phylogroups with ClermonTyping ####

"""
In this script we execute the script 'clermonTyping.sh' to identify E. coli phylogroups.
We use as input nucleotide fasta files (.fna), which must be given to the
ClermonTyping script with parallelization.
"""

import subprocess
from pathlib import Path
import pandas as pd

#=========================================================================
# Step 0: Stablish paths and input files
#=========================================================================
# Directories
genome_dir = Path("escherichia_ncbi_dataset").resolve()
results_dir = Path("results_clermon").resolve()
results_dir.mkdir(parents=True, exist_ok=True)

# ClermonTyping script
clermon_script = Path("./ClermonTyping-master/clermonTyping.sh").resolve()

#=========================================================================
# Step 1: Run clermonTyping
#=========================================================================
clermon_cmd = f"""
find "{genome_dir}" -type f -name "*.fna" -print0 | \
parallel -0 -j 32 --bar '
genome={{/.}};
mkdir -p "{results_dir}/$genome";
cd "{results_dir}/$genome";
"{clermon_script}" --fasta {{}} --name "$genome";
find . -type f ! -name "*phylogroups.txt" -delete
'
"""
subprocess.run(clermon_cmd, shell=True, check=True)

#=========================================================================
# Step 2: Create combined text file
#=========================================================================
# Read phylogroups files
phylogroups_files = sorted(results_dir.rglob("*phylogroups.txt"))

# Read dataframes
dfs = []
for file in phylogroups_files:
    try:
        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = [
            "fna_file",
            "genes",
            "genes_presence",
            "alleles",
            "phylogroup",
            "mash_group"
        ]
        df["assembly"] = "_".join(file.parent.name.split("_")[0:2])
        dfs.append(df)

    except Exception as e:
        print(f"Skipping {file}: {e}")

if not dfs:
    raise Exception("No phylogroups file found")

# Combine and export
combined_df = pd.concat(dfs, ignore_index=True, sort=False)
combined_df.to_csv(f"{results_dir}/all_phylogroups.txt", sep="\t", index=False)
print(f"Phylogroups saved to: {results_dir}/all_phylogroups.txt")

#=========================================================================
# Step 3: A little more cleanup
#=========================================================================
# 1. Count subdirectories in results directory (1 per genome)
subdirs = [d for d in results_dir.iterdir() if d.is_dir()]
n_subdirs = len(subdirs)

# 2. Count rows in combined file
n_rows = combined_df.shape[0]

print(f"Subdirectories found: {n_subdirs}")
print(f"Total rows: {n_rows}")

# 3. Conditional cleanup
if n_rows == n_subdirs:
    print("Congrats.")

    for item in results_dir.iterdir():
        # Remember to keep some output
        if item.name == "all_phylogroups.txt":
            continue
        # Remove subdirectories
        if item.is_dir():
            import shutil
            shutil.rmtree(item)
        # Just in case
        elif item.is_file():
            item.unlink()
    print(f"{n_subdirs} subdirectories removed.")
    print("Clean as a whistle. Only all_phylogroups.txt remains.")

else:
    print("Too bad!! The number of rows does not match the number of subdirectories.")
    print("Some genomes missing. All files remain.")