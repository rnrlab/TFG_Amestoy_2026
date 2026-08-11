##### Parsing HHblits results #####

"""
Parsing of results functional annotation with HHblits
over alignments of clusters' proteins, constructed
for a Pangenome with PanTA.
"""

import glob
import re

#=========================================================
# Step 1 - INPUTS
#=========================================================
base_dir = "hhblits_results"

hhr_files = glob.glob(f"{base_dir}/*.hhr")

#=========================================================
# Step 2 - PARSE
#=========================================================
results = {}

for hhr_file in hhr_files:
    cluster_name = hhr_file.split("/")[-1].split(".")[0]

    with open(hhr_file, "r") as f:
        lines = f.readlines()

    selected_hit = None

    for i, line in enumerate(lines):
        if line.startswith(">"):
            header = line.strip()
            next_line = lines[i + 1] if i + 1 < len(lines) else ""

            # Parse hit number
            hit_number = None
            for j in range(i-1, max(i-10, 0), -1):
                if lines[j].startswith("No "):
                    hit_number = lines[j].split()[1]
                    break

            # Parse header
            hit_clean = header[1:]
            parts = [p.strip() for p in hit_clean.split(";")]

            accession = parts[0]
            name = parts[1]
            description = ";".join(parts[2:])

            # Parse E-value
            evalue_hh = None
            match = re.search(r"E-value=([\deE\.\-]+)", next_line)
            if match:
                evalue_hh = match.group(1)

            # Filter hit
            bad_keywords = ["unknown", "uncharacterized", "hypothetical"]

            is_duf = "DUF" in accession or "DUF" in name
            is_unknown = any(k in description.lower() for k in bad_keywords)

            if not is_duf and not is_unknown:
                selected_hit = {
                    "hit_number": hit_number,
                    "accession": accession,
                    "name": name,
                    "description": description,
                    "evalue_hh": evalue_hh
                }
                break

    # si no encontró ninguno bueno
    if selected_hit is None:
        selected_hit = {
            "hit_number": None,
            "accession": None,
            "name": None,
            "description": None,
            "evalue_hh": None
        }

    results[cluster_name] = selected_hit

#=========================================================
# Step 3 - OUTPUTS
#=========================================================
out_tsv = "out_pipolin_clusters_annotations.tsv"
with open(out_tsv, "w") as f:
    f.write("gene\thit_number\taccession\tname\tdescription\tevalue_hh\n")
    for cluster_name, h in results.items():
        desc = h['description'].replace("\t", " ").replace("\n", " ")

        f.write(
            f"{cluster_name}\t"
            f"{h['hit_number']}\t"
            f"{h['accession']}\t"
            f"{h['name']}\t"
            f"{desc}\t"
            f"{h['evalue_hh']}\n"
        )


