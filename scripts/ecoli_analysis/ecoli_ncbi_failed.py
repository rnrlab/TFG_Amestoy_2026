#### Search for assemblies of failed genomes during NCBI-datasets donwload ####

"""
This script comes as consecuence of ecoli_ncbi_download, since some of the genomes
gave an error while downloading the selected files (fastas and GFF).

In this script we will look for those genome's assemblies, with the objective of
re-running the datasets download for them only.
"""

import json
import os

#=========================================================================
#===== Step 0: Reading input JSONL file ==================================
#=========================================================================

base_dir = "./escherichia_ncbi_dataset/ncbi_dataset/data"

all_assemblies = []

with open(f"{base_dir}/assembly_data_report.jsonl", encoding="utf-8") as in_file:
    for line in in_file:
        assembly = json.loads(line)
        all_assemblies.append(assembly['currentAccession'])

#=========================================================================
#===== Step 1: Export all assemblies =====================================
#=========================================================================
with open("ecoli_assemblies.txt", "w", encoding="utf-8") as out_file:
    for acc in all_assemblies:
        out_file.write(acc + "\n")

print(f'Total genomes: {len(all_assemblies)}')

#=========================================================================
#===== Step 2: Export failed assemblies ==================================
#=========================================================================
failed_assemblies = []
not_protein_assemblies = []

sub_dirs = set(os.listdir(base_dir))

for acc in all_assemblies:
    acc_path = f"{base_dir}/{acc}"

    if acc not in sub_dirs:
        failed_assemblies.append(acc)
    elif len(os.listdir(acc_path)) != 3:
        not_protein_assemblies.append(acc)

with open("ecoli_failed_assemblies.txt", "w", encoding="utf-8") as out_file_2:
    for acc in failed_assemblies:
        out_file_2.write(acc + "\n")
print(f'Failed genomes: {len(failed_assemblies)}')

with open("ecoli_not_protein_assemblies.txt", "w", encoding="utf-8") as out_file_3:
    for acc in not_protein_assemblies:
        out_file_3.write(acc + "\n")
print(f'Genomes w/o proteins: {len(not_protein_assemblies)}')