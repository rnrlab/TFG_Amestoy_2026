#### Screening of pipolins in bacterial genomes ####

"""
This script uses a hmm piPolB profile to look for piPolBs (that is, pipolins) genomes (in our case in bacterial genomes).
For this purpose, we first use the command-line tool NCBI-datasets to download dehydrated genomes (protein files)
of our selected taxon (in this case E. coli --> taxon_id: 2).

This dehydration is usually necessary for data storage purposes.
Batch processing is also common when processing more than 200 orders to the NCBI.

Results of the HMMSearch are saved in a single tabulated text file.
"""

import subprocess
import os
import re

#=========================================================================
# Step 0: Extract the dehydrated genomes and decompress the directory
#=========================================================================
subprocess.run("./datasets download genome taxon 2 --assembly-source GenBank --assembly-version latest --exclude-atypical --include protein --dehydrated", shell=True)
subprocess.run("unzip ncbi_dataset.zip", shell=True)

#=========================================================================
# Step 1: Save original dehydrated genomes to another file
#=========================================================================
with open("ncbi_dataset/fetch.txt", "r") as genome_dh_file:
    genome_dh = genome_dh_file.readlines()
with open("ncbi_dataset/fetch_original.txt", "w") as genome:
    genome.write("".join(genome_dh))

#=========================================================================
# Step 2: Run HMMSearch and save results
#=========================================================================
with open("PiPol_Coinc.txt", "w") as pipol_coinc:
    # Iterate to rehydrate the proteomes in groups of 100
    for n in range(0,len(genome_dh),100):
        genomes_proteins = {}
        group = genome_dh[n:n+100]
        # Rewrite fetch.txt with the group of 100
        with open("ncbi_dataset/fetch.txt", "w") as fetch_group:
            fetch_group.write("".join(group))
        # Rehydrate the proteomes, save them in a single .faa file, and perform hmmsearch
        subprocess.run("./datasets rehydrate --max-workers 30 --directory ./", shell=True)
        for root, dirs, files in os.walk("ncbi_dataset/data"):
            for file in files:
                if file.endswith(".faa"):
                    dir_name = os.path.basename(root)
                    with open(os.path.join(root, file), "r") as protein_file:
                        protein_file_r = protein_file.read()
                    protein_name = re.compile(r'^>(\S+)', re.MULTILINE)
                    proteins = re.findall(protein_name, protein_file_r)
                    for protein in proteins:
                        genomes_proteins[protein] = dir_name
        subprocess.run("cat ncbi_dataset/data/*/protein.faa > ncbi_dataset/data/all_protein.faa", shell=True)
        subprocess.run("hmmsearch --tblout result.tbl --domtblout result.domtblout -E 1e-50 --cpu 30 pipolb_infomap_modified.hmm ncbi_dataset/data/all_protein.faa", shell=True)
        # Write output into results tabulated text file
        with open("result.domtblout", "r") as result_file:
            for line in result_file:
                if "#" not in line:
                    line_list = line.split()
                    data1 = "\t".join(line_list[0:22])
                    data2 = "_".join(line_list[22:])
                    protein = line.split()[0]
                    assembly = genomes_proteins[protein]
                    pipol_coinc.write(assembly + "\t" + data1 + "\t" + data2 + "\n")
        # Remove subdirectories with proteins (storage purposes)
        subprocess.run('rm -r "ncbi_dataset/data"/GCA*', shell=True)

#=========================================================================
# Step 3: Add header and save
#=========================================================================
subprocess.run("cat encabezado_corto.txt PiPol_Coinc.txt > results_hmmsearch.txt", shell=True)



