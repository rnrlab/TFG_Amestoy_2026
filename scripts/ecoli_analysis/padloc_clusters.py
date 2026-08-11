#### Identification of Defense Systems genes within protein clusters ####

"""
As inputs we use the multifasta files of proteins within each cluster,
selected from a Pangenome with PanTA.

We compare these files with the hmm profiles from PADLOC-db.
"""

import os
import subprocess
import pandas as pd
from tqdm import tqdm
import time

clusters_dir = "out_pipolin_clusters_fastas"
hmms_dir = "padloc-db-master/hmm"
outdir = "results_hmm_padloc"
"""
os.makedirs(outdir, exist_ok=True)

commands = []

#==========================================================
# Step 1 - Generate commands
#==========================================================
clusters = [f for f in os.listdir(clusters_dir) if f.endswith(".fasta")]
hmms = [h for h in os.listdir(hmms_dir) if h.endswith(".hmm")]

for fasta in tqdm(clusters, desc="Generando comandos"):
    fasta_path = os.path.join(clusters_dir, fasta)
    base = fasta.replace(".fasta", "")
    cluster_out = os.path.join(outdir, base)
    os.makedirs(cluster_out, exist_ok=True)

    for hmm in hmms:
        hmm_path = os.path.join(hmms_dir, hmm)
        hmm_base = os.path.basename(hmm)

        tblout = os.path.join(cluster_out, f"{hmm_base}.tbl")
        logout = os.path.join(cluster_out, f"{hmm_base}.log")

        cmd = (
            f"hmmsearch --cpu 4 "
            f"--tblout {tblout} "
            f"{hmm_path} {fasta_path} > {logout}"
        )

        commands.append(cmd)

total_jobs = len(commands)
print(f"Total jobs: {total_jobs}")

#==========================================================
# Step 2 - Execute GNU Parallel
#==========================================================
proc = subprocess.Popen(
    ["parallel", "-j", "8", "--joblog", "parallel.log"],
    stdin=subprocess.PIPE,
    text=True
)

proc.stdin.write("\n".join(commands))
proc.stdin.close()

#==========================================================
# Step 3 - Progress bar
#==========================================================
pbar = tqdm(total=total_jobs, desc="Ejecutando hmmsearch")

done = 0
while proc.poll() is None:  # mientras parallel siga vivo
    if os.path.exists("parallel.log"):
        with open("parallel.log") as f:
            done = sum(1 for line in f) - 1  # restar cabecera
        pbar.n = done
        pbar.refresh()
    time.sleep(1)

pbar.close()
print("Done executing.")
"""
#==========================================================
# Step 4 - Parsing results
#==========================================================
def parse_hmm_metadata(hmm_file):
    meta = {
        "hmm_name": None,
        "hmm_accession": None,
        "hmm_description": None,
        "hmm_length": None,
        "system": None,
    }

    with open(hmm_file) as f:
        for line in f:
            if line.startswith("NAME"):
                meta["hmm_name"] = line.split()[1]
            elif line.startswith("ACC"):
                meta["hmm_accession"] = line.split()[1]
            elif line.startswith("DESC"):
                meta["hmm_description"] = " ".join(line.split()[1:])
            elif line.startswith("LENG"):
                meta["hmm_length"] = int(line.split()[1])

    # Infer system from NAME
    name = meta["hmm_name"]
    if name:
        prefix = name.split("_")[0].lower()
        meta["system"] = prefix

    return meta

rows = []

for cluster in os.listdir(outdir):
    cluster_path = os.path.join(outdir, cluster)
    if not os.path.isdir(cluster_path):
        continue

    for tbl in os.listdir(cluster_path):
        if not tbl.endswith(".tbl"):
            continue

        hmm = tbl.replace(".tbl", "")
        tbl_path = os.path.join(cluster_path, tbl)
        hmm_path = os.path.join(hmms_dir, hmm)

        # Parse metadata from the HMM file
        meta = parse_hmm_metadata(hmm_path)

        with open(tbl_path, "r") as f:
            for line in f:
                if line.startswith("#") or line.strip() == "":
                    continue

                parts = line.split()
                if len(parts) < 8:
                    continue

                target = parts[0]
                evalue = float(parts[4])
                score = float(parts[5])
                bias = float(parts[6])

                rows.append({
                    "cluster": cluster,
                    "hmm": hmm,
                    "target": target,
                    "evalue": evalue,
                    "score": score,
                    "bias": bias,

                    # Add HMM metadata
                    "hmm_name": meta["hmm_name"],
                    "hmm_accession": meta["hmm_accession"],
                    "hmm_description": meta["hmm_description"],
                    "hmm_length": meta["hmm_length"],
                    "system": meta["system"],
                })


#==========================================================
# Step 5 - Complete DataFrame
#==========================================================
df = pd.DataFrame(rows)
df.to_csv(f"{outdir}/padloc_all_hits.csv", index=False)
print("Total hits:", len(df))

#==========================================================
# Step 6 - Best hit per cluster
#==========================================================
df_best = (df
           .sort_values(["cluster", "evalue", "score"], ascending=[True, True, False])
           .groupby("cluster")
           .head(1)
           .reset_index(drop=True))

df_best = df_best.rename(columns={
    "hmm": "best_hmm",
    "target": "best_target",
    "evalue": "best_evalue",
    "score": "best_score",
})

df_best.to_csv(f"{outdir}/padloc_best_hits_per_cluster.csv", index=False)

#==========================================================
# Step 7 - Presence-absence matrix
#==========================================================
presence = (df
            .groupby(["cluster", "hmm"])
            .size()
            .unstack(fill_value=0)
            .astype(bool)
            .astype(int))

presence = presence.sort_index()
presence = presence.reindex(sorted(presence.columns), axis=1)

presence.to_csv(f"{outdir}/padloc_hmm_matrix.csv")
