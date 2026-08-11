#### Parse nucleotide fasta files for PPanGGOLiN ####

"""
To work with nucleotide fasta files (.fna) as input, PPanGGOLiN requires a TSV file
with a specific structure:
    column_1    column_2
    <assembly>  <path_to_fasta_file>
With 1 row per genome.
"""

import os

#=========================================================================
#===== Step 0: Stablish paths ============================================
#=========================================================================
base_dir = "./escherichia_ncbi_dataset"
out_dir = "./ecoli_ppanggolin_input.fasta.list"

#=========================================================================
#===== Step 1: Save assemblies and fasta files paths =====================
#=========================================================================
assemblies_paths = {}
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".fna"):
            assembly = os.path.basename(root)
            file_path = os.path.abspath(os.path.join(root, file))
            assemblies_paths[assembly] = file_path

#=========================================================================
#===== Step 2: Write TSV file ============================================
#=========================================================================
with open(out_dir, "w") as outfile:
    for assembly in assemblies_paths:
        outfile.write(f"{assembly}\t{assemblies_paths[assembly]}\n")