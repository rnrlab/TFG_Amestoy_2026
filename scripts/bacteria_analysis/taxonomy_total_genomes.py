#### Taxonomy Annotation from an Assembly for the Complete Bacteria Genomes ####

"""
Script very similar to taxonomy_from_assembly
Except that in this script we will start from a pre-existing JSON file: assembly_data_report.jsonl
It is automatically created when screening with datasets (see HMMSearch_Final)
"""

import ujson
from tqdm import tqdm

#=========================================================================
#===== Step 0: Stablish paths and define variables =======================
#=========================================================================

input_file = "ncbi_dataset/data/assembly_data_report.jsonl"
output_file = "total_genomes_taxonomy.txt"

# To save memory, we'll write to the output_file every 10000 lines
output_lines = [] # buffer: list to store the results in groups of 10000
count = 0
buffer_size = 10000 # Buffer size

#=========================================================================
#===== Step 1: Extractint metadata from json =============================
#=========================================================================

with open(input_file, "r", encoding="utf-8") as f:
    with open(output_file, "w") as out:
        # Header
        out.write("assembly\torganism_name\tspecies\tgenus\n")

        for line in tqdm(f, desc="Processing assemblies"): # tqdm para ver el progreso
            # Load data and genome accession
            data = ujson.loads(line)
            accession = data["currentAccession"]

            # Extract organism name and clean it
            organism_name = data['organism']['organismName']
            organism_name_clean = organism_name.replace("[", "").replace("]", "").replace("'", "").replace(" ", "_")
            words = organism_name_clean.split("_")

            # Extract genus
            genus = words[0]

            # If the organism name only has the genus (1 word), we consider species = genus
            if len(words) >= 2:
                species = "_".join(words[:2])
            else:
                species = genus

            # Save results in a list
            output_lines.append(f"{accession}\t{organism_name_clean}\t{species}\t{genus}\n")
            count += 1

            # Write results in block for each buffer_size
            if count % buffer_size == 0:
                out.writelines(output_lines)
                output_lines = []

        # We write any remaining lines into the buffer
        if output_lines:
            out.writelines(output_lines)

print(f"Procesadas {count} líneas. Archivo guardado en {output_file}.")
