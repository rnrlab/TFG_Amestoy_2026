#### Extracting Information from ExplorePipolin Results ####

"""
We use the ExplorePipolin results saved in a directory as input.
This directory contains a subdirectory for each genome.

This script is a more complete version of info_from_EP_v2.
"""

import os
import glob
from Bio import SeqIO

# ecoli: "ecoli_results_EP"
# Bacteria: "EP_results_bact"
base_dir = "EP_results_bact"
pipolin_dictionary = {}  # format: pipolin_id : [values]

for genome_folder in glob.glob(f"{base_dir}/*"):
    genome_id = genome_folder.replace(f"{base_dir}/", "")

    if not os.path.exists(f"{base_dir}/{genome_id}/{genome_id}/pipolins"):
        print(f"There is not pipolins folder for {genome_id}")

        pipolin_id = genome_id.split(".")[0]

        pipolin_dictionary[pipolin_id] = [
            genome_id,  # assembly
            genome_id,  # Genome_ID
            '0',        # pipolin_presence
            'NA',       # Pipolin_file
            'NA'        # Pipolin_scaffold
        ]

        remaining_columns = 25 #!!!! --> will need to change this number if columns change
        pipolin_dictionary[pipolin_id] += ['NA'] * remaining_columns

        continue

    for pipolin in os.listdir(f"{base_dir}/{genome_id}/{genome_id}/pipolins"):
        if "single_record" not in pipolin and "v0" in pipolin and ".gbk" in pipolin:
            pipolin_id = pipolin.split(".")[0].replace("_v0", "")

            # Record [pipolin_id](key-Field 0), [genome_id], [genome_num], [pipolin_presence], [pipolin file], [reconstruction type]
            pipolin_dictionary[pipolin_id] = [genome_id]
            pipolin_dictionary[pipolin_id].append(genome_id.replace("G_", ""))
            pipolin_dictionary[pipolin_id].append('1')
            pipolin_dictionary[pipolin_id].append(pipolin)
            pipolin_dictionary[pipolin_id].append(pipolin.split(".")[1])

            # Record  [number of piPolBs], [list of piPolB lengths], [sum of piPolB lengths], [pipolin length multirecord], [pipolin fragments]
            num_pipolb = 0
            length_pipolb_list = []
            length_pipolb_sum = 0
            pipolin_length = 0 # Suma directa de las longitudes de los fragmentos de pipolina (.gbk)
            pipolin_fragments = 0
            piPolB_check_fragmented = False
            for record in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin}", "genbank"):
                pipolin_fragments += 1
                pipolin_length += len(record.seq)

                for pipolin_cds in record.features:
                    if pipolin_cds.type == "CDS":
                        if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                            num_pipolb += 1
                            length_pipolb_list.append(len(pipolin_cds.qualifiers["translation"][0]))
                            length_pipolb_sum += len(pipolin_cds.qualifiers["translation"][0])
                            piPolB_check_fragmented = True

            pipolin_dictionary[pipolin_id].append(str(num_pipolb))
            pipolin_dictionary[pipolin_id].append(str(length_pipolb_list))
            pipolin_dictionary[pipolin_id].append(str(length_pipolb_sum))
            pipolin_dictionary[pipolin_id].append(str(pipolin_length))
            pipolin_dictionary[pipolin_id].append(str(pipolin_fragments))

            #Record [pipolin length in single record] and [length check]
            pipolin_length_reconstruction = 0 # Suma reconstruida de las longitudes de los fragmentos de la pipolina (.single_record.gbk)
            piPolB_check_reconstructed = False
            for record_rec in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin.replace(".gbk", ".single_record.gbk")}", "genbank"):
                pipolin_length_reconstruction = len(record_rec.seq)
                pipolin_cds_list_rec = [feat_rec for feat_rec in record_rec.features if feat_rec.type == "CDS"]
                for pipolin_cds_rec in pipolin_cds_list_rec:
                    if pipolin_cds_rec.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                        piPolB_check_reconstructed = True
            pipolin_dictionary[pipolin_id].append(str(pipolin_length_reconstruction)) # Anota la suma reconstruida de las longitudes
            pipolin_dictionary[pipolin_id].append(str(pipolin_length+100*(pipolin_fragments-1) == pipolin_length_reconstruction)) # Booleano = True, siempre que la suma directa sea aprox. igual a la suma reconstruida

            #check if [piPolBs in fragmented?], [pipolb in reconstruction?] [10]
            pipolin_dictionary[pipolin_id].append(str(piPolB_check_fragmented))
            pipolin_dictionary[pipolin_id].append(str(piPolB_check_reconstructed))

            #get pipolin contig identifiers
            pipolin_fa = pipolin.replace(".gbk", ".fa")
            pipolin_contig_ids = [re.id for re in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin_fa}", "fasta")]

            #Record pipolin fragments length [11], pipolin fragments length sum [12], pipolin contigs IDs [13], pipolin coordinates in contig [14]
            pipolin_fragment_length = [] # Longitud de los fragmentos de la pipolina
            pipolin_fragment_length_sum = 0 # Esto debería ser igual a pipolin_length
            pipolin_contig = [] # IDs de los contigs --> MUY IMPORTANTE
            pipolin_contig_coordinates = [] # Coordenadas de la(s) pipolina(s) EN EL CONTIG: ">contig_accesion inicio:final"
            for fa in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin_fa}", "fasta"):
                if fa.id in pipolin_contig_ids:
                    pipolin_fragment_length.append(str(len(fa)))
                    pipolin_fragment_length_sum += len(fa)
                    pipolin_contig.append(str(fa.id))
                    pipolin_contig_coordinates.append(fa.description.split(" ")[1]) # ">contig_accesion inicio:final"

            pipolin_dictionary[pipolin_id].append(str(pipolin_fragment_length))
            pipolin_dictionary[pipolin_id].append(str(pipolin_fragment_length_sum))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_coordinates))

            #Record [repeats presence], [att number], [coordinates], [att lengths], [att mean length], [att type], [att overlaps], [att overlaps clean]
            att_presence = False
            att_number = 0
            att_coordinates = []
            att_length = []
            att_mean_length = 0
            att_type = []
            att_overlaps = []
            att_overlaps_clean = []

            for record in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin}", "genbank"):
                repeats = [feat for feat in record.features if (feat.type == "repeat_region" and feat.qualifiers["rpt_family"][0]=="Att")]
                if repeats != []:
                    att_presence = True
                    att_number += len(repeats)

                    for repeat in repeats:
                        att_coordinates.append(str(record.id)+"::"+str(repeat.location))
                        att_length.append(int(abs(repeat.location.end-repeat.location.start)))

                        if repeat.qualifiers["note"][0] not in att_type:
                            att_type.append(repeat.qualifiers["note"][0])

                        att_range = range(int(repeat.location.start), int(repeat.location.end))
                        for feat in record.features:
                            if feat.type != "repeat_region" and feat.type != "source":
                                if feat.location.start in att_range or feat.location.end in att_range:
                                    if "product" in list(feat.qualifiers.keys()):
                                        att_overlaps.append(
                                            str(record.id) + "::" + str(feat.qualifiers["product"]).replace("'", ""))
                                        if str(feat.qualifiers["product"]).replace("'", "").replace("[", "").replace(
                                                "]", "") not in att_overlaps_clean:
                                            att_overlaps_clean.append(
                                                str(feat.qualifiers["product"]).replace("'", "").replace("[",
                                                                                                         "").replace(
                                                    "]", ""))
                                    else:  # In case it has no product
                                        att_overlaps.append(str(record.id) + "::" + str(feat.type))
                                        if str(feat.type) not in att_overlaps_clean:
                                            print(str(feat.type), att_overlaps_clean)
                                            att_overlaps_clean.append(
                                                str(feat.type).replace("'", "").replace("[", "").replace("]", ""))

            pipolin_dictionary[pipolin_id].append(str(att_presence))
            pipolin_dictionary[pipolin_id].append(str(att_number))
            pipolin_dictionary[pipolin_id].append(str(att_type))
            pipolin_dictionary[pipolin_id].append(str(att_coordinates))
            pipolin_dictionary[pipolin_id].append(str(att_length))
            if len(att_length) > 0:
                pipolin_dictionary[pipolin_id].append(str(sum(att_length)/len(att_length)))
            else:
                pipolin_dictionary[pipolin_id].append("0")
            pipolin_dictionary[pipolin_id].append(str(att_overlaps))
            if att_overlaps_clean == []:
                pipolin_dictionary[pipolin_id].append("NA")
            else:
                pipolin_dictionary[pipolin_id].append(", ".join(att_overlaps_clean))


            #Record [assembly gaps multirecord], [assembly gaps single_record], [assembly gaps_paired-ends], [assembly gaps_pipolin_Structure]
            reconstruction_gaps_multirecord = 0
            for record in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin}", "genbank"):
                for pipolin_feat in record.features:
                    if "assembly_gap" in str(pipolin_feat.type):
                        reconstruction_gaps_multirecord += 1
            pipolin_dictionary[pipolin_id].append(str(reconstruction_gaps_multirecord))

            reconstruction_gaps_single_record = 0
            gap_paired_ends = 0
            gap_pipolin_structure = 0
            for record in SeqIO.parse(f"{base_dir}/{genome_id}/{genome_id}/pipolins/{pipolin.replace(".gbk", ".single_record.gbk")}", "genbank"):
                for pipolin_feat in record.features:
                    if "assembly_gap" in str(pipolin_feat.type):
                        reconstruction_gaps_single_record += 1

                        if pipolin_feat.qualifiers["linkage_evidence"][0] == "paired-ends":
                            gap_paired_ends += 1

                        if pipolin_feat.qualifiers["linkage_evidence"][0] == "pipolin_structure":
                            gap_pipolin_structure += 1
            pipolin_dictionary[pipolin_id].append(str(reconstruction_gaps_multirecord == gap_paired_ends))
            pipolin_dictionary[pipolin_id].append(str(reconstruction_gaps_single_record))
            pipolin_dictionary[pipolin_id].append(str(gap_paired_ends))
            pipolin_dictionary[pipolin_id].append(str(gap_pipolin_structure))

#add haeder
table_output = "Pipolin_ID\tassembly\tGenome_ID\tpipolin_presence\tPipolin_file\tPipolin_scaffold\tpiPolB_num\tpiPolB_fragments_length\tpiPolB_fragments_length_sum\t"
table_output += "Pipolin_length\tPipolin_fragments\tPipolin_length_reconstruction\tLength_check\tpiPolB_check_fragments\tpiPolB_check_reconstruction\t"
table_output += "Pipolin_fragment_lengths\tPipolin_fragment_lengths_sum\tPipolin_contig_IDs\tPipolin_contig_coordinates\tAtt_presence\t"
table_output += "Att_number\tAtt_type\tAtt_coordinates\tAtt_lengths\tAtt_mean_length\tIntegration_site\tIntegration_site_clean\tAssembly_gaps_multirecord\t"
table_output += "Same_assembly_gaps_multi_and_paired_ends\tAssembly_gaps_single_record\tAssembly_gaps_paired_ends\tAssembly_gaps_pipolin_structure\n"

for key in pipolin_dictionary:
    table_output += str(key)+"\t"+"\t".join(pipolin_dictionary[key])+"\n"

with open("ecoli_info_from_EP.txt", "w") as f:
    f.write(table_output)
