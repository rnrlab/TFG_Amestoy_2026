##### Obtaining multifasta of proteins within clusters form PanTA #####

"""
This script looks for the protein FASTA sequences of each protein
within a cluster, then merges them into a multi-FASTA file.
"""

from Bio import SeqIO
import json
import os
from collections import defaultdict
from multiprocessing import Pool

#=========================================================
# Step 1 - INPUTS
#=========================================================
with open("results_panta/20_core_clusters.txt") as f:
    out_clusters = [line.strip() for line in f]

clusters_json = json.load(open("results_panta/annotated_clusters.json"))

multifasta_dir = "20_core_clusters"
os.makedirs(multifasta_dir, exist_ok=True)

#=========================================================
# Step 2 - FUNCIÓN PARA UN CLUSTER
#=========================================================
def process_cluster(c):

    records_to_write = []

    genes = clusters_json[c]['gene_id']

    assemblies = defaultdict(set)

    for g in genes:
        assembly = g.split("-")[0]
        assemblies[assembly].add(g)

    for assembly, gene_set in assemblies.items():

        fasta_path = f"results_panta/samples/{assembly}/{assembly}.faa"

        if not os.path.exists(fasta_path):
            continue

        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in gene_set:
                records_to_write.append(record)

    output_path = f"{multifasta_dir}/{c}.fasta"

    with open(output_path, "w") as out_f:
        SeqIO.write(records_to_write, out_f, "fasta")

    return c

#=========================================================
# Step 3 - PARALELIZACIÓN
#=========================================================
if __name__ == "__main__":

    n_cores = 16

    with Pool(n_cores) as p:
        results = p.map(process_cluster, out_clusters)

    print("Done.")