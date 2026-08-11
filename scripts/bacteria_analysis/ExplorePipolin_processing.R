
#=========================================================================
#===== Step 0: Importing and Parsing data ================================
#=========================================================================
### A) Import
## A.1) Genomes w/ pipolins
pipolin_metadata_filtered4 <- read.table(
  file = "pipolin_metadata_filtered4.txt",
  sep = '\t',
  header = TRUE
)

## A.2) Info from ExplorePipolin: info_from_EP_v2.txt <-- corrected version from info_from_EP.txt
# First incomplete version
info_from_EP <- read.csv(
  file = "info_from_EP.txt",
  sep = "\t"
)

# Eliminate old EP columns from pipolin_metadata_filtered4
old_columns <- colnames(info_from_EP)[-1] # Quitando la columna de assembly
keep_columns <- setdiff(colnames(pipolin_metadata_filtered4), old_columns)
pipolin_metadata_filtered5 <- pipolin_metadata_filtered4[, keep_columns]

### B) Export
# Export pipolin_metadata_filtered5 and merge in Python (faster)
write.table(
  pipolin_metadata_filtered5,
  file = "pipolin_metadata_filtered5.txt",
  sep = "\t",
  row.names = FALSE
)
