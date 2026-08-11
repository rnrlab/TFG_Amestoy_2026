##### Merging core genome clusters #####

"""
This script takes alignments of proteins of a given number
of core genome clusters (from a pangenome from PanTA)
and merges them all into one single file,
with the maximum possible number of different genomes.

The base directory only contains .aln.fasta files
"""

import os
from Bio import SeqIO

#=========================================================
# Step 0 - INPUTS
#=========================================================

multifasta_dir = "20_core_clusters_mafft"

cluster_files = sorted(os.listdir(multifasta_dir))

#=========================================================
# Step 1 - OBTAIN GENOMES
#=========================================================

all_genomes = set()

for cluster_file in cluster_files:
    file_path = os.path.join(multifasta_dir, cluster_file)

    for record in SeqIO.parse(file_path, "fasta"):
        genome = record.id.split("-")[0]
        all_genomes.add(genome)

print(f"Total genomes: {len(all_genomes)}")

#=========================================================
# Step 2 - PARSE
#=========================================================

all_genomes = sorted(all_genomes)
genomes_dic = {genome: [] for genome in all_genomes}

for cluster_file in cluster_files:
    file_path = os.path.join(multifasta_dir, cluster_file)

    ## A) Cluster dictionary
    # Take each genome only 1 time per cluster
    seen_genomes = set()

    # Temporary dictionary for cluster
    cluster_sequences = {}

    for record in SeqIO.parse(file_path, "fasta"):
        genome = record.id.split("-")[0]

        # Only 1 protein per genome
        if genome in seen_genomes:
            continue

        seen_genomes.add(genome)
        cluster_sequences[genome] = str(record.seq)

    ## B) Global dictionary
    # Concatenate sequences
    if len(cluster_sequences) == 0:
        continue

    length = len(next(iter(cluster_sequences.values())))

    for genome in genomes_dic:
        if genome in cluster_sequences:
            genomes_dic[genome].append(cluster_sequences[genome])

        # Add gaps if genome not present in this cluster
        else:
            genomes_dic[genome].append("-" * length)

#=========================================================
# Step 3 - CONCATENATED MULTI-FASTA
#=========================================================

out_file = "ecoli_concatenated_proteins.aln.fasta"

with open(out_file, "w") as out:
    for genome, seq_list in genomes_dic.items():
        full_seq = "".join(seq_list)
        out.write(f">{genome}\n{full_seq}\n")

print("Done.\n")