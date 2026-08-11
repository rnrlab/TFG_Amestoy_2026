#### Initial Processing of HMMSearch Results ####

"
Here we import the screening results table. 
Then we process it, removing duplicates and other errors.
"

#=========================================================================
#===== Step 0: Importing data from HMMSearch =============================
#=========================================================================

## A) Import
results_hmmsearch <- read.table("../results_hmmsearch.txt", header = TRUE, sep = "\t")

## B) Analysis and filtering of duplicates
# Number of genomes w/ and w/o duplicates
length(results_hmmsearch$assembly)

(duplicated_assemblies <- results_hmmsearch$assembly[duplicated(results_hmmsearch$assembly)])
length(duplicated_assemblies)

length(unique(results_hmmsearch$assembly))

# Number of proteins w/ and w/o duplicates
length(results_hmmsearch$target_name)

(duplicated_prots <- results_hmmsearch$target_name[duplicated(results_hmmsearch$target_name)])
length(duplicated_prots)

length(unique(results_hmmsearch$target_name))

# Remove different columns in duplicates
results_hmmsearch_filtered <- results_hmmsearch[, -c(11:23)]

# Remove duplicated proteins
results_hmmsearch_filtered <- unique(results_hmmsearch_filtered)

# Check duplicates absence
length(unique(results_hmmsearch_filtered$target_name))
length(unique(results_hmmsearch_filtered$assembly))

#=========================================================================
#===== Step 1: Exporting data from HMMSearch =============================
#=========================================================================

## A) Exporting filtered table
write.table(
  results_hmmsearch_filtered,
  file = "results_hmmsearch_filtered.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## B) Exporting only genomes assemblies
write.table(
  unique(results_hmmsearch_filtered$assembly),
  file = "pipol_genomes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#=========================================================================
#===== Step 2: Importing data of all bacterial genomes ===================
#=========================================================================

## A) Importing all-bacterial-genomes metadata table
gene_counts_summary <- read.table(
  "../gene_counts_summary.txt",
  header = TRUE,
  sep = "\t"
)

## B) Process coding and non-coding genes
# Genomes without coding genes (protein_coding == NA)
length(gene_counts_summary$assembly[is.na(gene_counts_summary$protein_coding)])

# Genomes with coding genes (protein_coding != NA)
length(gene_counts_summary$assembly[!is.na(gene_counts_summary$protein_coding)])

# We removed genomes without coding genes
gene_counts_summary <- gene_counts_summary[!is.na(gene_counts_summary$protein_coding), ]

#=========================================================================
#===== Step 3: Exporting data of all bacterial genomes ===================
#=========================================================================

## A) Exporting filtered table
write.table(
  gene_counts_summary,
  file = "gene_counts_summary_filtered.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## B) Exporting only genomes assemblies
write.table(
  unique(gene_counts_summary$assembly),
  file = "total_genomes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

