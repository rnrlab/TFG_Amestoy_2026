#### Analysis of clusters position relative to the pipolin ####

"
We work with the metadata from PanTA, from 2 pangenomes:
 - One for E. coli genomes
 - Another for E.coli pipolin-containing genomes
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
pipolins_clusters <- read.table(
  file = "panta_pipolins_clusters/t_gene_presence_absence_filtered_pipolins.Rtab",
  sep = "\t",
  header = T
)
colnames(pipolins_clusters)[colnames(pipolins_clusters) == "X"] <- "assembly"

panta_clusters_filtered_pipolins <- read.table(
  file = "panta_pipolins_clusters/panta_clusters_filtered_pipolins.txt",
  sep = "\t",
  header = T
)

ecoli_dataset_phylogroup <- read.table(
  file = 'ecoli_dataset_phylogroups.txt',
  sep = '\t',
  header = T
)
pipolins_clusters$pipolb <- ecoli_dataset_phylogroup$pipolin_presence[match(pipolins_clusters$assembly, ecoli_dataset_phylogroup$assembly)]

# Results from BLASTp alignment
blast <- read.table(
  file ="panta_pipolins_clusters/blast_best_per_subject.tsv",
  sep = "\t",
  header = T
)
colnames(blast)[colnames(blast) == "evalue"] <- "evalue_blast"

#=========================================================================
# Step 1 - Parse input data
#=========================================================================
## 1.1) General merge and analysis 
# Merge BLASTp results: EColi-Pangenome X Pipol-Pangenome
colnames(signif_df)[colnames(signif_df) == "evalue"] <- "evalue_interpro"

signif_df <- merge(signif_df, blast,
           by.x = "gene",
           by.y = "sseqid",
           all.x = TRUE)

# Analyze merged-dataset
length(unique(signif_df$qseqid))
sum(signif_df$evalue_blast <= 1e-20, na.rm = TRUE)
sum(signif_df$evalue_blast > 1e-20, na.rm = TRUE)

## 1.2) Filter table
# Eliminate redundant columns --> always keep the original
signif_df <- signif_df %>%
  select(
    -qcovs,
    -qcovhsp,
    -mismatch,
    -gapopen,
    -qstart,
    -qend,
    -sstart,
    -send
  )

# Calculate Long-Coverage and Short-Coverage
signif_df <- signif_df %>%
  mutate(
    qcov = (length / qlen) * 100,
    scov = (length / slen) * 100
  ) %>%
  relocate(qcov, scov, .before = evalue_blast)

# Filter by E-vale, Ident%, Query-Coverage and Subject-Coverage
signif_in_pipolin <- signif_df %>%
  filter(
    if_all(c(evalue_blast, pident, qcov, scov), ~!is.na(.)),
    evalue_blast <= 1e-20,
    pident >= 70,
    qcov >= 50,
    scov >= 50
  )

signif_out_pipolin <- signif_df %>%
  filter(
    if_any(c(evalue_blast, pident, qcov, scov), is.na) |
      evalue_blast > 1e-20 |
      pident < 70 |
      qcov < 50 |
      scov < 50
  )

#=========================================================================
# Step X - Export data
#=========================================================================
#write.table(signif_df, file = "panta_clusters/ecoli_significant_clusters_dataset.txt", sep = "\t", quote = F, row.names = F)
#write.table(signif_in_pipolin, file = "panta_clusters/ecoli_signif_in_pipolin_dataset.txt", sep = "\t", quote = F, row.names = F)
#write.table(signif_out_pipolin, file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", quote = F, row.names = F)
#writeLines(signif_out_pipolin$gene, con = "panta_pipolins_clusters/out_pipolin_clusters.txt")


