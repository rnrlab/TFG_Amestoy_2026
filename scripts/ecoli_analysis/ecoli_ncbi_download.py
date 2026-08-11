#### Searching for FASTA sequences of genomes and proteins and GFF files ####

"""
We use datasets only for Escherichia coli (taxon id: 562).
We download the dehydrated genomes and rehydrate them in groups.
VERY similar to HMMSearch_Final.py, but with other options for the datasets command
"""

import subprocess

#=========================================================================
#===== Step 0: Download dehydrated genomes with NCBI-datasets ============
#=========================================================================

# Include 'genome' (nucleotides fasta) and 'protein' (aas fasta)
subprocess.run("datasets download genome taxon 562 --include genome,protein,gff3 --assembly-source GenBank --assembly-version latest --exclude-atypical --assembly-level scaffold,chromosome,complete --dehydrated --filename escherichia_ncbi_dataset.zip", shell=True)
subprocess.run("unzip escherichia_ncbi_dataset.zip -d escherichia_ncbi_dataset", shell=True)

#=========================================================================
#===== Step 1: Save dehydrated genomes ===================================
#=========================================================================

# Move all the dehydrated genomes to a new permanent file (fetch_original.txt)
with open("escherichia_ncbi_dataset/ncbi_dataset/fetch.txt", "r") as dh_file:
    dh_genomes = dh_file.readlines()
with open("escherichia_ncbi_dataset/ncbi_dataset/fetch_original.txt", "w") as genomes:
    genomes.write("".join(dh_genomes))

#=========================================================================
#===== Step 2: Rehydrate genomes with error logging ======================
#=========================================================================
failed_genomes = []

# Iterate in groups of 100, rehydrating the genomes
for n in range(0,len(dh_genomes),100):
    group = dh_genomes[n:n+100]
    with open("escherichia_ncbi_dataset/ncbi_dataset/fetch.txt", "w") as fetch_group:
        fetch_group.write("".join(group))

    # "rehydrate" actúa directamente sobre fetch.txt --> cada fasta se guarda en un directorio cuyo nombre es el Assembly del genoma
    try:
        subprocess.run("datasets rehydrate --max-workers 30 --directory escherichia_ncbi_dataset", shell=True)
    except Exception as e:
        print(f"Batch {n} failed: {e}. Saving accessions to failed_genomes.txt")
        failed_genomes.extend(group)

#=========================================================================
#===== Step 3: Save failed genomes to another file =======================
#=========================================================================
if failed_genomes:
    with open("escherichia_ncbi_dataset/ncbi_dataset/failed_fetch.txt", "w") as f:
        f.writelines(failed_genomes)
        print(f"{len(failed_genomes)} failed. Retry with failed_genomes.txt")
else:
    print("All genomes rehydrated successfully!!!")
