#### ExplorePipolin Data Processing ####

"
Here we work with output files of (almost) all previous R scripts.
The point is to process data in order to make some pretty plots later,
in another script though.
"

#=========================================================================
#===== Step 0: Importing data ============================================
#=========================================================================

## A) Table of prevalence analysis by genus
bacteria_summary <- read.table(
  file = "bacteria_summary.txt",
  header = TRUE,
  sep = "\t"
)

## B) Table of assemblages with their organism - species - genus
pipol_genomes_taxonomy <- read.table(
  file = "pipol_genomes_taxonomy.txt",
  header = TRUE,
  sep = "\t"
)

## C) CheckM2 Completeness and Contamination Table
checkm2 <- read.table(
  file = "checkm2_all_genomes.tsv",
  sep = "\t",
  header = TRUE
)

## D) Table with metadata of EP results
info_from_EP <- read.table(
  file = "info_from_EP.txt",
  sep = "\t",
  header = TRUE
)

#=========================================================================
#===== Step 1: Tables management =========================================
#=========================================================================

## A) Take assembly only
checkm2$assembly <- sapply(strsplit(checkm2$Name, "_"), function(x) paste(x[1:2], collapse="_"))

## B) Change columns names
colnames(bacteria_summary)[1:4] <- c("genus", "family", "order", "class")

## C) Merge tables to take genomes assemblies and their taxonomy
# 1st merge
pipolin_metadata <- merge(pipol_genomes_taxonomy, bacteria_summary[1:4], by = "genus")

# 2nd merge
pipolin_metadata <- pipolin_metadata[, c("assembly", "organism_name", "species", "genus", "family", "order", "class")]

# 3rd merge - WRONG -
pipolin_metadata_2 <- merge(pipolin_metadata, checkm2, by = "assembly")

# Difference of 1 assembly (16 417 --> 16 416)
setdiff(pipolin_metadata$assembly, pipolin_metadata_2$assembly)
pipolin_metadata$assembly[grep("GCA_044987585", pipolin_metadata$assembly)]
checkm2$assembly[grep("GCA_044987585", checkm2$assembly)]

# Wrong version: GCA_044987585.1 --> GCA_044987585.2
pipolin_metadata$assembly[pipolin_metadata$assembly == "GCA_044987585.1"] <- "GCA_044987585.2"

# 3rd merge - CORRECTED -
pipolin_metadata <- merge(pipolin_metadata, checkm2, by = "assembly")
rm(pipolin_metadata_2) # This was just to prove what was wrong

# 4th (and last) merge
pipolin_metadata <- merge(pipolin_metadata, info_from_EP, by = "assembly")

## D) Processing genomes w/o pipolins
genomes_without_pipolins <- pipolin_metadata[pipolin_metadata$num_of_pipolins == 0,]

# Check genomes
cat('Genomes w/o pipolins:', nrow(genomes_without_pipolins))

# Check genera
cat('Genus with genomes w/o pipolins:', length(unique(genomes_without_pipolins$genus)))
genomes_without_pipolins_genus <- as.data.frame(table(genomes_without_pipolins$genus))

# Check Escherichia and Vibrio --> great mayority so it is no problem
cat('Escherichia genomes w/o pipolins (%):', 
    (genomes_without_pipolins_genus$Freq[genomes_without_pipolins_genus$Var1=="Escherichia"]/sum(genomes_without_pipolins_genus$Freq))*100
)
cat('Vibrio genomes w/o pipolins (%):', 
    (genomes_without_pipolins_genus$Freq[genomes_without_pipolins_genus$Var1=="Vibrio"]/sum(genomes_without_pipolins_genus$Freq))*100
)
cat('Escherichia and Vibrio genomes w/o pipolins (%):', 
    (genomes_without_pipolins_genus$Freq[genomes_without_pipolins_genus$Var1=="Escherichia"]/sum(genomes_without_pipolins_genus$Freq) +
      genomes_without_pipolins_genus$Freq[genomes_without_pipolins_genus$Var1=="Vibrio"]/sum(genomes_without_pipolins_genus$Freq))*100
)

# Filter pipolin_metadata of genomes without pipolins
genomes_with_pipolins <- pipolin_metadata[pipolin_metadata$num_of_pipolins != 0,]

#=========================================================================
#===== Step 2: Exporting data ============================================
#=========================================================================

# Genomes w/ pipolins detected by HMMSearch (16 417)
write.table(
  pipolin_metadata,
  file = "pipolin_metadata.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Genomes w/ pipolins detected by HMMSearch and ExplorePipolin (15 987)
write.table(
  genomes_with_pipolins,
  file = "pipolin_metadata_filtered.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
