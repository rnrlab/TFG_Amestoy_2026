#### Identification of Escherichia coli sequence types (ST) with Multi-locus sequence typing (MLST) ####

"""
This script executes mlst to identify sequence types or STs between E. coli genomes.
For this, the input will be a FASTA file (.fna) for each genome.
As output, mlst in --full mode returns a TSV file with the following columns:
 - FILE: input filename
 - SCHEME: auto-detected scheme
 - ST: sequence type assigned
 - STATUS: quality of genotype
 - SCORE: score of genotype
 - ALLELES: identified alleles
"""

from pathlib import Path
import subprocess
import pandas as pd

#=========================================================================
# Step 0: Stablish base and output directories
#=========================================================================
# Directories
genome_dir = Path('PRUEBA_escherichia_ncbi_dataset')
results_dir = Path('results_mlst')
results_dir.mkdir(parents=True, exist_ok=True)

#=========================================================================
# Step 1: Execute MLST
#=========================================================================
mlst_cmd = f'''
find "{genome_dir}" -type f -name "*.fna" -print0 | \
parallel -0 -j 32 --bar '
fna="{{}}"
assembly=$(basename "$(dirname "$fna")")
outfile="{results_dir}/$assembly.txt"
mlst --full --label "$fna" --outfile "$outfile" "$fna"
'
'''
subprocess.run(mlst_cmd, shell=True, check=True)

#=========================================================================
# Step 2: Create combined text file
#=========================================================================
# Read all output text files
all_outputs = sorted(results_dir.glob("*.txt"))

# Read tables
dfs = []
for file in all_outputs:
    try:
        df = pd.read_csv(file, sep="\t")
        df['assembly'] = file.stem
        dfs.append(df)
    except Exception as e:
        print(f"Skipping {file}: {e}")

if not dfs:
    raise Exception("No output files found")

# Combine and export
combined_df = pd.concat(dfs, ignore_index=True)
combined_df.to_csv(results_dir / "results_mlst.txt", sep="\t", index=False)
print(f"Results written to {results_dir}/results_mlst.txt")

#=========================================================================
# Step 3: Some cleanup
#=========================================================================
# 1. Count singel text files (1 per genome)
text_files = list(results_dir.glob("GCA*.txt"))
n_text_files = len(text_files)

# 2. Count rows in combined file
n_rows = combined_df.shape[0]

print(f"Single text files found: {n_text_files}")
print(f"Total rows: {n_rows}")

# 3. Conditional cleanup
if n_text_files == n_rows:
    print("Congrats.")
    for file in text_files:
        file.unlink()
    print(f"{n_text_files} single text files removed.")
    print("Clean as a whistle. Only results_mlst.txt remains.")

else:
    print("Too bad!! The number of rows does not match the number of single text files.")
    print("Some genomes missing. All files remain.")