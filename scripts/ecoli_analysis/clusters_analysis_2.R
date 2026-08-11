#### Analysis of clusters annotation from InterProScan ####

"
We work with the metadata resulting of runnning InterProScan
with the FAA of the selected clusters from PanTA.
"

#=========================================================================
# Step 0 - Import and parse input data
#=========================================================================
### A) Import
## Presence-absence
clusters <- read.table(
  file = "panta_clusters/t_gene_presence_absence_filtered.Rtab",
  sep = "\t",
  header = T
)
colnames(clusters)[colnames(clusters) == "X"] <- "assembly"

ecoli_dataset_phylogroup <- read.table(
  file = "ecoli_dataset_phylogroups.txt",
  sep = "\t",
  header = T
)

# Adding manual piPolB cluster
clusters$pipolb <- ecoli_dataset_phylogroup$pipolin_presence[match(clusters$assembly, ecoli_dataset_phylogroup$assembly)]

## InterProScan
signif_df <- read.table(
  file = "panta_clusters/ecoli_significant_clusters_dataset.txt",
  sep = "\t",
  header = T
)

significant_clusters_best_hit <- read.csv(
  file = "panta_clusters/ecoli_significant_clusters_best_hit.txt",
  sep = "\t",
  header = T,
  row.names = NULL
)

### B) Parse
# Missing clusters --> not annotated by InterProScan
setdiff(signif_df$gene, significant_clusters_best_hit$gene)

# Merge dataframes
signif_df <- merge(signif_df, significant_clusters_best_hit, by = "gene")

#=========================================================================
# Step 1 - Positive / Negative correlation analysis
#=========================================================================
# Create contingency tables and add to dataset
contingency_ls <- lapply(clusters[, !(names(clusters) %in% c("assembly", "pipolb"))], function(col) {
  t <- table(clusters$pipolb, col)
  data.frame(
    pNcN = t[1,1],
    pNcY = t[1,2],
    pYcN = t[2,1],
    pYcY = t[2,2]
  )
})

contingency_df <- do.call(rbind, contingency_ls)
contingency_df$gene <- rownames(contingency_df)
rownames(contingency_df) <- NULL

signif_df <- merge(signif_df, contingency_df, by = "gene")

#=========================================================================
# Step X - Export data
#=========================================================================
#write.table(signif_df, file = "panta_clusters/ecoli_significant_clusters_dataset.txt", sep = "\t", quote = F, row.names = F)
#write.table(signif_df$gene, file = "panta_clusters/ecoli_significant_clusters_ID.txt", sep = "\t", quote = F, row.names = F, col.names = F)
