#### Obtaining the taxonomy of genomes with and without pipolins ####

"
We combined the results of HMMSearch (results_hmmsearch_filtered.txt)
with a file containing the organism, species, and genus of the assemblies (pipol_genomes_taxonomy.txt). 
The HMMSearch results have duplicate assemblies, while the taxonomy results do not.
We did the same for complete bacterial genomes.
"

packages <- c("dplyr", "taxize")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
#===== Step 0: Importing data ============================================
#=========================================================================

## A) Genomes w/ pipolins info and initial taxonomy
results_hmmsearch_filtered <- read.table(
  file = "results_hmmsearch_filtered.txt", 
  header = TRUE, 
  sep = "\t")

pipol_genomes_taxonomy <- read.table(
  file = "pipol_genomes_taxonomy.txt",
  header = TRUE,
  sep = "\t"
)

## B) Total genomes (w/ and w/o pipolins) info and initial taxonomy
gene_counts_summary <- read.table(
  file = "gene_counts_summary_filtered.txt",
  header = TRUE,
  sep = "\t"
)

total_genomes_taxonomy <- read.delim(
  file = "total_genomes_taxonomy.txt",
  header = TRUE,
  sep = "\t"
)

# We only take genomes w/ proteins
total_genomes_taxonomy <- total_genomes_taxonomy[
  total_genomes_taxonomy$assembly %in% gene_counts_summary$assembly,
]

#=========================================================================
#===== Step 1: Frequency of genomes and Prevalence =======================
#=========================================================================

## A) Frequency of genomes w/ pipolins in each genus
pipol_summary <- as.data.frame(sort(table(pipol_genomes_taxonomy$genus), decreasing = TRUE))
colnames(pipol_summary) <- c("Genus", "Pipol.Genomes")
pipol_summary$Pipol.Prop <- pipol_summary$Pipol.Genomes/sum(pipol_summary$Pipol.Genomes)

## B) Frequency of total genomes in each genus
total_summary <- as.data.frame(sort(table(total_genomes_taxonomy$genus), decreasing = TRUE))
colnames(total_summary) <- c("Genus", "Total.Genomes")

## C) Processing
# Merge tables
bacteria_summary <- merge(pipol_summary, total_summary, by = "Genus")

# Proportion: Genomes.with.pipolins / Total.genomes
bacteria_summary$Genomes.Prop <- bacteria_summary$Pipol.Genomes/bacteria_summary$Total.Genomes

# Prevalence by genus
bacteria_summary$Pipol.Prevalence <- bacteria_summary$Genomes.Prop*100

#=========================================================================
#===== Step 2: Obtaining taxonomy (taxize library) =======================
#=========================================================================

## A) Extract Genera and search taxonomy
genera <- unique(as.character(bacteria_summary$Genus))

taxonomy <- classification(genera, db = "ncbi")

## B) Extract Order
orders <- sapply(taxonomy, function(x) {
  if(is.data.frame(x)) {
    order <- x$name[x$rank == "order"]
    if(length(order) == 0) return(NA)
    return(order)
  } else {
    return(NA)
  }
})

# Change to dataframe
df_order <- data.frame(Genus = names(orders), Order = orders, stringsAsFactors = FALSE)

# Add to summary table
bacteria_summary <- left_join(bacteria_summary, df_order, by="Genus")
bacteria_summary <- bacteria_summary %>%
  relocate(Order, .after = Genus)

## C) Extract Family
families <- sapply(taxonomy, function(x) {
  if(is.data.frame(x)) {
    family <- x$name[x$rank == "family"]
    if(length(family) == 0) return(NA)
    return(family)
  } else {
    return(NA)
  }
})

# Change to dataframe
df_family <- data.frame(Genus = names(families), Family = families, stringsAsFactors = FALSE)

# Add to summary table
bacteria_summary <- left_join(bacteria_summary, df_family, by="Genus")
bacteria_summary <- bacteria_summary %>%
  relocate(Family, .after = Genus)

## D) Extract Class
classes <- sapply(taxonomy, function(x) {
  if(is.data.frame(x)) {
    class <- x$name[x$rank == "class"]
    if(length(class) == 0) return(NA)
    return(class)
  } else {
    return(NA)
  }
})

# Change to dataframe
df_class <- data.frame(Genus = names(classes), Class = classes, stringsAsFactors = FALSE)

# Add to summary table
bacteria_summary <- left_join(bacteria_summary, df_class, by="Genus")
bacteria_summary <- bacteria_summary %>%
  relocate(Class, .after = Order)

#=========================================================================
#===== Step 3: Output processing =========================================
#=========================================================================

## A) Processing NAs
bacteria_summary$Family[is.na(bacteria_summary$Family)] <- "Unknown"
bacteria_summary$Order[is.na(bacteria_summary$Order)] <- "Unknown"
bacteria_summary$Class[is.na(bacteria_summary$Class)] <- "Unknown"

## B) Export summar table
write.table(
  bacteria_summary, 
  file = "bacteria_summary.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)
