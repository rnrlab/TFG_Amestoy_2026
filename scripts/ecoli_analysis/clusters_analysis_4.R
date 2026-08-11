#### Analysis of clusters annotation by HHblits ####

"
We work with the results of functional annotation of clusters
with HHblits.

Refering only to the clsuters defined (by us) ot be outside 
the pipolin.
"

packages <- c("dplyr")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
# Step 0 - Import input data
#=========================================================================
# Dataset of clusters outside the pipolin --> includes annotation w/ InterProScan
signif_out_pipolin <- read.table(
  file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt",
  sep = "\t",
  header = T
)

# Results from HHblits
hhblits <- read.table(
  file = "panta_clusters/out_pipolin_clusters_annotations.tsv",
  sep = "\t",
  header = T,
)

#=========================================================================
# Step 1 - Parse input data
#=========================================================================
# Process information from InterProScan
colnames(signif_out_pipolin)[
  colnames(signif_out_pipolin) %in% c(
    "best_hit", "accession", "description")
] <- c("name_interpro", "accession_interpro", "description_interpro")

# Process information from HHblits
colnames(hhblits)[
  colnames(hhblits) %in% c(
    "accession", "name", "description")
] <- c("accession_hh", "name_hh", "description_hh")

# Merge datasets
signif_out_pipolin <- merge(signif_out_pipolin, hhblits, by = "gene")

#=========================================================================
# Step 2 - Export
#=========================================================================
write.table(signif_out_pipolin, file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", quote = F, row.names = F)
