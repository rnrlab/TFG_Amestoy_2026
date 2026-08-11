#### Searching for GFF files starting from the GenBank Assembly (GCA_XXXXXXX) ####

"""
We use datasets starting from a text file with genome assemblies --> pipol_genomes.txt --> this could easily be changed to choose the file name.
We download the dehydrated genomes and rehydrate them in groups.
VERY similar to HMMSearch_Final
"""

import subprocess

#=========================================================================
#===== Step 0: Download dehydrated genomes with NCBI-datasets ============
#=========================================================================

# Include 'genome' (nucleotides fasta) and 'protein' (aas fasta)
subprocess.run("datasets download genome accession --inputfile pipol_genomes.txt --include gff3 --assembly-source GenBank --assembly-version latest --dehydrated --filename gff_ncbi_dataset.zip", shell=True)
subprocess.run("unzip gff_ncbi_dataset.zip -d gff_ncbi_dataset", shell=True)

#=========================================================================
#===== Step 1: Save dehydrated genomes ===================================
#=========================================================================

# Move all the dehydrated genomes to a new permanent file (fetch_original.txt)
with open("gff_ncbi_dataset/ncbi_dataset/fetch.txt", "r") as dh_file:
    dh_genomes = dh_file.readlines()
with open("gff_ncbi_dataset/ncbi_dataset/fetch_original.txt", "w") as genomes:
    genomes.write("".join(dh_genomes))

#=========================================================================
#===== Step 2: Rehydrate genomes =========================================
#=========================================================================

# Iterate in groups of 100, rehydrating the genomes
for n in range(0,len(dh_genomes),100):
    group = dh_genomes[n:n+100]
    with open("gff_ncbi_dataset/ncbi_dataset/fetch.txt", "w") as fetch_group:
        fetch_group.write("".join(group))
    # "rehydrate" actúa directamente sobre fetch.txt --> cada fasta se guarda en un directorio cuyo nombre es el Assembly del genoma
    subprocess.run("datasets rehydrate --max-workers 30 --directory gff_ncbi_dataset", shell=True)
