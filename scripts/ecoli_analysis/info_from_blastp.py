#### Extracting information from BLASTp results ####

"""

"""


import pandas as pd

# Cargar datos
cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qcovs", "qcovhsp", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

blast_tsv = pd.read_csv(
    "results_panta_pipolin/blast_clusters_results.tsv",
    sep="\t",
    names=cols
)


# Ordenar por calidad del hit
blast_sorted = blast_tsv.sort_values(
    by=["sseqid", "evalue", "bitscore"],
    ascending=[True, True, False]
)

# Quedarse con el mejor hit por cada subject
blast_best = blast_sorted.drop_duplicates(subset="sseqid", keep="first")

print(len(blast_best))  # debería ser 266

# Guardar resultado
blast_best.to_csv(
    "results_panta_pipolin/blast_best_per_subject.tsv",
    sep="\t",
    index=False
)
