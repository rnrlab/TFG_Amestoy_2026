#### Extracting information from BLASTp results ####

"""
This script extracts and organizes results from a BLASTp search,
as a method of functional annotation of the protein clusters.
"""


import pandas as pd

# Load data
cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qcovs", "qcovhsp", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

blast_tsv = pd.read_csv(
    "results_panta_pipolin/blast_clusters_results.tsv",
    sep="\t",
    names=cols
)


# Order by quality of hit
blast_sorted = blast_tsv.sort_values(
    by=["sseqid", "evalue", "bitscore"],
    ascending=[True, True, False]
)

# Keep only the best hit per subject
blast_best = blast_sorted.drop_duplicates(subset="sseqid", keep="first")

print(len(blast_best)) # just checkin'

# Save results
blast_best.to_csv(
    "results_panta_pipolin/blast_best_per_subject.tsv",
    sep="\t",
    index=False
)
