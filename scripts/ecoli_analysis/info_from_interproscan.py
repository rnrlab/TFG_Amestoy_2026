##### Parse output information from InterProScan #####

"""
We use as input de JSON file provided by InterProScan as output
(just one of the outputs of said software).
"""

import json

#=========================================================================
# Step 1 - Write JSON output
#=========================================================================
# Read input data
data = json.load(open('panta_clusters/ecoli_significant_clusters_part3.json', 'r'))

# Write dictionary and write summary
matches_for_cluster = {}
for protein in data['results']:
    cluster_complete_name = protein['xref'][0]['name']
    cluster_name = protein['xref'][0]['id']
    print(f'{cluster_complete_name}:')

    lines = []
    matches_for_cluster[cluster_name] = {}

    for i, m in enumerate(protein['matches']):
        try:
            match = m['signature']['entry']['name']
            accession = m['signature']['entry']['accession']
            description = m['signature']['entry']['description']
        except TypeError:
            match = m['signature']['name']
            accession = m['signature']['accession']
            description = m['signature']['description']

        try:
            evalue = m['evalue']
        except KeyError:
            evalue = None

        matches_for_cluster[cluster_name][match] = {}
        matches_for_cluster[cluster_name][match]['accession'] = accession
        matches_for_cluster[cluster_name][match]['description'] = description
        matches_for_cluster[cluster_name][match]['evalue'] = evalue

        line = f" {i}) {match} | {accession} | {description} | {evalue} |"
        lines.append(line)

        if len(lines) == 0:
            lines.append("Not a single match")

    print("\n".join(lines))
    print("\n")

# Write output json
with open('panta_clusters/info_from_interproscan_part3.json', 'w') as f:
    json.dump(matches_for_cluster, f, indent=4)

#=========================================================================
# Step 2 - Write TSV output
#=========================================================================
# Filter proteins by e-value
data2 = json.load(open('panta_clusters/info_from_interproscan_part3.json', 'r'))

genes = {}
for gene in data2:
    options = [
        (match_name, match_info['accession'], match_info['description'], match_info['evalue'])
        for match_name, match_info in data2[gene].items()
    ]
    genes[gene] = options

best_hits = {}
for gene, matches in genes.items():
    valid = [
        (name, accession, description, e)
        for name, accession, description, e in matches
        if e is not None
    ]

    if valid:
        best_name, best_accession, best_description, best_evalue = min(valid, key=lambda x: x[3])
    else:
        if matches:
            best_name, best_accession, best_description, best_evalue = matches[0]
        else:
            best_name = best_accession = best_description = best_evalue = None

    best_hits[gene] = {
        "best_hit": best_name,
        "accession": best_accession,
        "description": best_description,
        "best_evalue": best_evalue
    }

# Write output TSV
with open('panta_clusters/ecoli_significant_clusters_part3_best_hit.txt', 'w') as f:
    f.write("gene\tbest_hit\taccession\tdescription\tevalue\n")
    for gene in best_hits:
        f.write(f"{gene}\t"
                f"{best_hits[gene]['best_hit']}\t"
                f"{best_hits[gene]['accession']}\t"
                f"{best_hits[gene]['description']}\t"
                f"{best_hits[gene]['best_evalue']}\n")
