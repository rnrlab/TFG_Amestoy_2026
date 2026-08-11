#### Extracting information from IntegronFinder2.0 results ####

"""
We use the results from integron_finder (IntegronFinder2.0) as input.

The base directory contains a subdirectory for each genome, which contains 3 main files:
    - mysequences.integrons: file with all integrons and their elements detected in all sequences in the input file.
    - mysequences.summary: summary file with the number and type of integrons per sequence.
    - integron_finder.out: copy standard output.
The base directory name goes as: "Results_Integron_Finder_mysequences",
where "mysequences" can be (e.g.): "GCA_000225665.2_ASM22566v2_genomic" (i.e. the genome assembly)

For each section, the input file is enclosed in square brackets [...].
"""

import os
import pandas as pd
import json
import logging

#=========================================================================
# Step 0: Define input directory
#=========================================================================

base_directory = "results_integronfinder"
out_text_file = "info_from_integronfinder.txt"
out_json_file = "info_from_integronfinder.json"
logging.basicConfig(
    filename="info_from_integronfinder.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

#=========================================================================
# Step 1: Extract data
#=========================================================================

### A) Number of genomes w/ and w/o integrons, and number of integrons per genome [mysequences.integrons, mysequences.summary]
G_list = []
G_w_info = []
G_w_integrons = {}
G_wo_integrons = []
other_files = []

# Look across base directory
for root, dirs, files in os.walk(base_directory):
    ## A.1) Genome assembly
    for n, directory in enumerate(dirs):
        assembly = "_".join(directory.split("_")[3:5])
        G_list.append(assembly)

        for file in os.listdir(root + '/' + directory):
            file_path = root + '/' + directory + '/' + file

            ## A.2) Number of genomes w/ and w/o integrons + output per genome [mysequences.integrons]
            if file.startswith("GCA") and file.endswith(".integrons"):
                G_w_info.append(assembly)

                # Separate genomes with and without integrons
                with open(os.path.join(file_path), "r") as integrons_file:
                    integrons_file = integrons_file.readlines()
                    if "# No Integron found" in integrons_file[2]: # See line 3 of "mysequences.integrons" to understand
                        logging.info(f'{n} {assembly} --> No integrons found :(')
                        G_wo_integrons.append(assembly)
                    else:
                        logging.info(f'{n} {assembly} --> INTEGRONS HERE!!! :)')
                        G_w_integrons[assembly] = {
                            "id_integron": {},
                            "total": 0,
                            "CALIN": 0,
                            "complete": 0,
                            "In0": 0
                        }

                        # Read integrons file
                        integrons_file_read = pd.read_csv(os.path.join(file_path), sep='\t', comment='#')

                        # Group by integrons
                        grouped_replicons = integrons_file_read.groupby('ID_integron')['ID_replicon'].unique().apply(list)

                        # Save info per integron in dictionary
                        for integron, replicon in grouped_replicons.items():
                            # Extract attCs info
                            attcs_df = integrons_file_read[
                                (integrons_file_read['ID_integron'] == integron) &
                                (integrons_file_read['element'].str.contains('attc', case=False))
                            ]

                            n_attcs = len(attcs_df)

                            attcs_dict = {}
                            for _, row in attcs_df.iterrows():
                                attcs_dict[row['element']] = {
                                    'pos_beg': row['pos_beg'],
                                    'pos_end': row['pos_end'],
                                    'distance_2attC': ('NA' if pd.isna(row['distance_2attC']) else row['distance_2attC']),
                                }

                            # Extract integron type
                            integron_type = integrons_file_read.loc[
                                integrons_file_read['ID_integron'] == integron, 'type'
                            ].iloc[0]

                            # Save to dictionary
                            G_w_integrons[assembly]["id_integron"][integron] = {
                                'type': integron_type,
                                'replicons': replicon,
                                'n_attC': n_attcs,
                                'attCs': attcs_dict
                            }

            else:
                other_files.append(file)

            ## A.3) Number of integrons per type and per genome [mysequences.summary]
            if file.startswith("GCA") and file.endswith(".summary"):
                if assembly in G_w_integrons:
                    # Read summary file
                    summary_file = pd.read_csv(os.path.join(file_path), sep='\t', comment='#')

                    num_calin = summary_file['CALIN'].sum()
                    num_complete_integron = summary_file['complete'].sum()
                    num_inO = summary_file['In0'].sum()

                    # Save info in dictionary
                    G_w_integrons[assembly]['total'] = int(num_calin) + int(num_complete_integron) + int(num_inO)
                    G_w_integrons[assembly]['CALIN'] = int(num_calin)
                    G_w_integrons[assembly]['complete'] = int(num_complete_integron)
                    G_w_integrons[assembly]['In0'] = int(num_inO)

#=========================================================================
# Step 2: Export data
#=========================================================================

# Export as text file
with open(out_text_file, "w") as text_file:
    text_file.write("assembly\t"
                    "total_integrons\t"
                    "calin\t"
                    "complete_integrons\t"
                    "In0\n")
    for G in G_w_integrons:
        text_file.write(f"{G}\t{G_w_integrons[G]['total']}\t{G_w_integrons[G]['CALIN']}\t{G_w_integrons[G]['complete']}\t{G_w_integrons[G]['In0']}\n")

# Export as json
with open(out_json_file, "w") as json_file:
    json.dump(G_w_integrons, json_file, indent = 2)

#=========================================================================
# Step 3: Summary to log file
#=========================================================================

logging.info("## COUNTS SUMMARY ##")
logging.info(f"Total genomes: {len(G_list)}")
logging.info(f"Genomes w/ info: {len(G_w_info)}")
logging.info(f" - Genomes w/o integrons: {len(G_wo_integrons)}")
logging.info(f" - Genomes w/ integrons: {len(G_w_integrons.keys())}")
logging.info(f"Other files: {len(other_files)}")
