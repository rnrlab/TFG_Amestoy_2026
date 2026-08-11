#### Analysis of ExplorePipolin results over 33k E. coli genomes ####

"
Our main input file is ecoli_info_from_EP.txt.
Take into account that most of these genomes will NOT have pipolins,
since we are parsing all of our E. coli genomes
"

packages <- c()
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
# Step 0: Import and parse input data
#=========================================================================
### A) Import
ecoli_info_from_EP <- read.csv(
  file = "ecoli_info_from_EP.txt",
  sep = "\t"
)

ecoli_pipolin_genomes <- readLines(
  con = "ecoli_pipolin_genomes.txt"
)

### B) Parse
# Count genomes w/ and w/o pipolins
nrow(ecoli_info_from_EP[ecoli_info_from_EP$pipolin_presence == '1',])
nrow(ecoli_info_from_EP[ecoli_info_from_EP$pipolin_presence == '0',])
nrow(ecoli_info_from_EP[ecoli_info_from_EP$pipolin_presence == '1',]) + nrow(ecoli_info_from_EP[ecoli_info_from_EP$pipolin_presence == '0',])

# Check new pipolin-containing-genomes
all(ecoli_pipolin_genomes %in% ecoli_info_from_EP$assembly)
length(setdiff(ecoli_info_from_EP$assembly[ecoli_info_from_EP$pipolin_presence == '1'], ecoli_pipolin_genomes))

#=========================================================================
# Step 1: Export data
#=========================================================================
# Save UPDATED ecoli genomes wiht pipolins
ecoli_pipolin_genomes_updated <- ecoli_info_from_EP$assembly[ecoli_info_from_EP$pipolin_presence == '1']
write.table(
  ecoli_pipolin_genomes_updated,
  file = 'ecoli_pipolin_genomes_updated.txt',
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

