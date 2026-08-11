#### Script to run Prodigal to create GFF3 files ####

"""
This script uses Prodigal to create a GFF3 file having an FNA file as input.
It starts from a directory with FNA files, one per genome, and runs
Prodigal in parallel.

This way prodigal creates a GFF3 file, then we manually add the "##FASTA" header
followed by the FNA sequence.
"""

import os
import subprocess

#=========================================================================
# Step 0: Stablish paths
#=========================================================================
base_dir = "escherichia_ncbi_dataset/ncbi_dataset/data"
dest_dir = "ecoli_gff_for_panta"
os.makedirs(dest_dir, exist_ok=True)

#=========================================================================
# Step 1: Run Prodigal over genomes
#=========================================================================
# With all genomes
prodigal_cmd = f"""
find {base_dir} -type f -name "*.fna" | \
parallel -j 32 --verbose '
assembly=$(basename "$(dirname "{{}}")")
prodigal -i {{}} -o {dest_dir}/$assembly.gff -f gff -p meta
'
"""
subprocess.run(prodigal_cmd, shell=True, check=True)

# In case the process is cut short [<file_with_missing_genomes>.txt] --> remove "#" from subprocess.run(...) below
prodigal_cmd_resume = f"""
cat ecoli_no_gff.txt | \
parallel -j 32 --verbose '
fna=$(echo {base_dir}/{{}}/*.fna)
prodigal -i "$fna" -o {dest_dir}/{{}}.gff -f gff -p meta
'
"""
#subprocess.run(prodigal_cmd_resume, shell=True, check=True)

#=========================================================================
# Step 3: Add FASTA to GFF3 files
#=========================================================================
for gff_file in os.listdir(dest_dir):
    # Get GFF
    gff_path = os.path.join(dest_dir, gff_file)
    genome_name = ".".join(gff_file.split('.')[0:2])

    # Search for FNA
    fna_path = None
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            if f.endswith(".fna") and genome_name in f:
                fna_path = os.path.join(root, f)
                break
        if fna_path:
            break

    # Copy ##FASTA to GFF
    if fna_path:
        with open(gff_path, "a") as gff_out, open(fna_path) as fna_in:
            # Copy ##FASTA header
            gff_out.write("##FASTA\n")

            # Copy FASAT sequence
            for line in fna_in:
                gff_out.write(line)
        print(f"Added FASTA to {gff_path}")
    else:
        print(f"No FNA found for {gff_file}")
