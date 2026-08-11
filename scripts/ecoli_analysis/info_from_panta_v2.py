#### Extracting information from PanTA results ####

"""
We use as input the results from PanTA, saved in a Rtab file.
This Rtab is a presence-absence matrix, with as many columns as genomes,
and as many rows as cluster groups.

We look to extract a summary text file with information about each found cluster,
as well as a text file with those cluster that are of interest to us.
"""

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

#=========================================================================
# Step 0: Define input file
#=========================================================================
rtab_file = 'results_panta_pipolin/gene_presence_absence.Rtab'
rtab_out = 'results_panta_pipolin/gene_presence_absence_filtered_pipolins.Rtab'

fasta_in = 'results_panta_pipolin/representative_clusters_prot.fasta'
fasta_out = 'results_panta_pipolin/representative_clusters_prot_filtered_pipolins.fasta'

#=========================================================================
# Step 1: Detect column names and prepare dtype
#=========================================================================
# Read columns
sample = pd.read_table(rtab_file, sep='\t', nrows=5)

# Separate clusters columns and convert the remaining in int16 (0|1)
cols = sample.columns
dtype_map = {col: "int16" for col in cols[1:]}

# Chunks
chunksize = 10000
summary_chunks = []

# Count lines for progress bar
num_lines = sum(1 for line in open(rtab_file))

#=========================================================================
# Step 2 - Parse file by chunks
#=========================================================================
for chunk in tqdm(
        pd.read_table(
            rtab_file,
            sep="\t",
            chunksize=chunksize,
            dtype=dtype_map,
            low_memory=False,
            engine="c"
        ),
        total=num_lines // chunksize,
        desc="Processing Rtab file"
    ):

    # Firs column: clusters names
    groups = chunk.iloc[:, 0]

    # Remaining columns (presence-absence)
    data = chunk.iloc[:, 1:]

    presence = (data > 0).sum(axis=1)
    absence = data.shape[1] - presence
    pct = (presence / data.shape[1] * 100).round(2)

    # CLassification --> IMPORTANTE REVISAR TAMBIEN !!!!!
    status = pd.cut(
        pct,
        bins=[-0.01, 1, 95, 99, 100],
        labels=["cloud", "shell", "soft_core", "core"],
        right=False
    )

    # Create data frame
    df = chunk.iloc[:, :1].copy()
    df["presence"] = presence
    df["absence"] = absence
    df["percentage"] = pct
    df["status"] = status

    summary_chunks.append(df)

#=========================================================================
# Step 3: Combine results
#=========================================================================
# Merge chunks
summary = pd.concat(summary_chunks, ignore_index=True)

# Save
summary.to_csv('results_panta_pipolin/panta_clusters_summary_pipolins.txt', sep="\t", index=False)

# Filter
filtered = summary[(summary["percentage"] >= 1) & (summary["percentage"] <= 100)] #Recordar cambiar si se quiere !!!!!!!!!
filtered.to_csv("results_panta_pipolin/panta_clusters_filtered_pipolins.txt", sep="\t", index=False)

print('Summary files saved.')
print("Done.")

#=========================================================================
# Step 4: Filter fasta file
#=========================================================================
print('Filtering FASTA now...')

# Load filtered clusters
filtered = pd.read_csv('results_panta_pipolin/panta_clusters_filtered_pipolins.txt', sep="\t")
groups_filter = set(filtered["Gene"])

# Filter fasta
records = (
    record for record in SeqIO.parse(fasta_in, "fasta")
    if record.id in groups_filter
)

# Output
count = SeqIO.write(records, fasta_out, "fasta")

print(f'Saved {count} records in {fasta_out}')
print("Done.")

#=========================================================================
# Step 5: Filter Rtab file
#=========================================================================
# 1) Read filtered groups
filtered = pd.read_csv("results_panta_pipolin/panta_clusters_filtered_pipolins.txt", sep="\t")
groups_filter = set(filtered["Gene"])   # MUY importante que sea un set

# 2) Separate clusters columns and convert the remaining in uint8 (0|1)
sample = pd.read_table(rtab_file, sep="\t", nrows=5)
cols = sample.columns
dtype_map = {col: "int16" for col in cols[1:]}

chunksize = 10000
num_lines = sum(1 for n in open(rtab_file))

# 3) Output file
with open(rtab_out, "w") as out:

    # Original header
    out.write("\t".join(cols) + "\n")

    # 4) Process by chunks
    for chunk in tqdm(
        pd.read_table(
            rtab_file,
            sep="\t",
            chunksize=chunksize,
            dtype=dtype_map,
            low_memory=False,
            engine="c"
        ),
        total=num_lines // chunksize,
        desc="Filtering Rtab"
    ):

        # Filter groups
        mask = chunk.iloc[:, 0].isin(groups_filter)
        filtered_chunk = chunk[mask].copy()

        # Binarize presence-absence -> now uint8 because
        filtered_chunk.iloc[:, 1:] = (filtered_chunk.iloc[:, 1:] > 0).astype("uint8")

        # Write output
        filtered_chunk.to_csv(out, sep="\t", index=False, header=False)

print("Filtered Rtab saved as:", rtab_out)
print("Done.")