#### Extracting information from CONJScan results ####

"""
We use the results from CONJScan, saved in a directory, as input.
This directory contains a subdirectory for each genome.

For each section, the input file is enclosed in square brackets [...].

CONJScan also detects MOB relaxases to identify plasmids, but it does so
including all MOB types in one category ("MOB"). Tipically you want to
differentiate between MOB types to separate plasmids types, but this
requires an aditional parsing (Step 1-B).
"""

import os
import pandas as pd

#=========================================================================
# Step 0: Define input directory
#=========================================================================
base_directory = "results_macsyfinder"
G_list = os.listdir(base_directory)

#=========================================================================
# Step 1: Extract data
#=========================================================================

### A) Extract number of systems detected in each genome [best_solution_summary.tsv]
all_genomes = []
for G in G_list:
    n_systems_path = os.path.join(base_directory, G, "best_solution_summary.tsv")
    n_systems_file = pd.read_csv(n_systems_path, sep='\t', comment='#')
    n_systems_file.insert(0, "assembly", G)
    all_genomes.append(n_systems_file)
systems_per_genome = pd.concat(all_genomes, ignore_index = True)
systems_per_genome.to_csv(f"{base_directory}/systems_per_genome.txt", sep = '\t', index = False)

### B) Extract MOB systems found per genome [best_solution.tsv]
G_to_MOB_systems = {}
for G in G_list:
    G_to_MOB_systems[G] = set()

    best_solution_path = os.path.join(base_directory, G, "best_solution.tsv")
    try:
        best_solution = pd.read_csv(best_solution_path, sep='\t', comment='#')
    except pd.errors.EmptyDataError:
        continue
    # Column 'gene_name' contains the HMM profile that detected the hit
    for value in best_solution['gene_name']:
        if 'MOB' in value:
            G_to_MOB_systems[G].add(value)
with open(f'{base_directory}/mobs_per_genome.txt', "w") as mobs_per_genome:
    mobs_per_genome.write("assembly\t"
                          "MOB_types\n")
    for G in G_list:
        mob_set = G_to_MOB_systems[G]

        mobs_per_genome.write(G + "\t")
        if not mob_set:
            mobs_per_genome.write("NA\n")
        else:
            mobs_per_genome.write(','.join(sorted(mob_set)) + "\n")

## B.1) Presence-Absence Matrix [mobs_per_genome.txt]
mobs_per_genome = pd.read_csv(f'{base_directory}/mobs_per_genome.txt', sep='\t')
mobs_matrix = mobs_per_genome['MOB_types'].str.get_dummies(sep=',')
mobs_matrix.insert(0, "assembly", mobs_per_genome['assembly'])
mobs_matrix.to_csv(f"{base_directory}/mobs_matrix.txt", sep = '\t', index = False)