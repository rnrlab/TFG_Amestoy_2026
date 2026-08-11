#### Taxonomy Annotation from an Assembly ####

"""
We will work with a text file containing genome assemblies.
We want to find the taxonomy of this genome, specifically: Organism, Species, Genus.

We will use datasets to download a JSON file with information about the genomes.
From this JSON file, we will extract the name of the organism to which it belongs.

We will finish by creating a summary table as output.
"""

import subprocess
import json

#=========================================================================
#===== Step 0: Obtaining json data =======================================
#=========================================================================

# Assembly input
input_assemblies = "pipol_genomes.txt"
output_json = "pipol_genomes_info.json"
output_taxonomy = "pipol_genomes_taxonomy.txt"

# Obtaining the json
command = f"./datasets summary genome accession --inputfile {input_assemblies} > {output_json}"
subprocess.run(command, shell=True)

#=========================================================================
#===== Step 1: Extractint metadata from json =============================
#=========================================================================

with open(output_json, "r", encoding="utf-8") as json_file:
    data = json.load(json_file)

    # Open results file (output in table format)
    with open(output_taxonomy, "w", encoding="utf-8") as taxonomy_file:
        taxonomy_file.write("assembly\torganism_name\tspecies\tgenus\n")

        # We iterate through the JSON, taking the organism name, species, and genus
        for report in data["reports"]:
            accession = report["accession"]
            organism_name = report["organism"]["organism_name"]
            organism_name_not_spaces = organism_name.replace(" ", "_")
            words = organism_name.split()

            genus = words[0]

            # If the organism name only has the genus (1 word), we consider species = genus
            if len(words) >= 2:
                species = "_".join(words[:2])
            else:
                species = genus
            taxonomy_file.write(f"{accession}\t{organism_name_not_spaces}\t{species}\t{genus}\n")
