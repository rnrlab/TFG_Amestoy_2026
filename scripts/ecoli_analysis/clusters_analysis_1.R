#### Analysis of clusters from PanTA over E. coli genomes ####

"
We work with the metadata of protein clusters found by PanTA over
32.700 E. coli genomes processed.

PanTA already has a cluster (and it may have several clusters actually) 
that identify with the piPolB protein. 
However, we do not care and will add our own cluster with the
pipolb presence-absence information we know from ExplorePipolin.
"

packages <- c("lsr", "dplyr", "psych")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
# Step 0 - Import and parse input data
#=========================================================================
### A) Import
clusters <- read.table(
  file = "panta_clusters/t_gene_presence_absence_filtered_original.Rtab",
  sep = "\t",
  header = T
)
colnames(clusters)[colnames(clusters) == "X"] <- "assembly"

ecoli_dataset_phylogroup <- read.table(
  file = 'ecoli_dataset_phylogroups.txt',
  sep = '\t',
  header = T
)

panta_clusters_filtered <- read.table(
  file = "panta_clusters/panta_clusters_filtered.txt",
  sep = "\t",
  header = T
)

setdiff(clusters$assembly, ecoli_dataset_phylogroup$assembly) # just checking, expecting "character(0)"

### B) Adding manual piPolB cluster
clusters$pipolb <- ecoli_dataset_phylogroup$pipolin_presence[match(clusters$assembly, ecoli_dataset_phylogroup$assembly)]
table(clusters$pipolb, clusters$groups_14801)

#=========================================================================
# Step 1 - Chi-Square and Cramer's V tests
#=========================================================================
### A) Run Chi-2 and CramersV tests. Store in list of lists.
results <- lapply(clusters[, !(names(clusters) %in% c("assembly", "pipolb"))], function(col) {
  tab <- table(clusters$pipolb, col)
  chi_test <- chisq.test(tab)
  cramer_test <- cramersV(tab)
  phi_test <- phi(tab)
  list(
    statistic = chi_test$statistic,
    p.value = chi_test$p.value,
    cramersV = cramer_test,
    phi = phi_test
  )
})

### B) Create dataframes from list of list. Bind dataframes.
results_df <- do.call(rbind, lapply(names(results), function(n) {
  data.frame(
    gene = n,
    chi.s = results[[n]]$statistic,
    p.value = results[[n]]$p.value,
    cramersV = results[[n]]$cramersV,
    phi = results[[n]]$phi
  )
}))
rownames(results_df) <- NULL

### C) Plot results
## Simple histogram --> range too big
hist(results_df$p.value,
     main = "Distribución de p-value",
     xlab = "p-value",
     col = "skyblue",
     border = "white",
     breaks = 30)

## -Log10 filtered histrogram
# Check distribution of Cramer's V
summary(results_df$cramersV)
min(results_df$p.value[results_df$p.value != 0]) # true minimum that R processes


# First we change 0s into the lowest p.value
p <- results_df$p.value
p[p == 0] <- min(results_df$p.value[results_df$p.value != 0])

# Now filter
vals <- -log10(p[results_df$cramersV > median(results_df$cramersV)]) # Change condition if wanted
length(vals)

hist(vals, 
     main = "Distribución de -log10(p-value)",
     xlab = "-log10(p-value)",
     col = "skyblue",
     border = "white")


#=========================================================================
# Step 2 - FDR correction
#=========================================================================
results_df$FDR <- p.adjust(results_df$p.value, method = "BH")

results_df <- results_df %>% relocate(FDR, .after = p.value)

#=========================================================================
# Step 3 - Filter data
#=========================================================================
### A) Genomes with piPolbs
sum(clusters$pipolb)
sum(clusters$pipolb)/nrow(clusters)*100

### B) Add number of genomes by cluster
results_df <- merge(results_df, panta_clusters_filtered, by.x = "gene", by.y = "Gene")

### C) Filter table
# First decide the filtering method
sum(results_df$FDR < 0.01) # NOT real filter here

# Moderate filtering
sum(results_df$cramersV > 0.15)
sum(results_df$cramersV > 0.15 & results_df$FDR < 0.01)

# Strong filtering
sum(results_df$cramersV > 0.20)
sum(results_df$cramersV > 0.20 & results_df$FDR < 0.01)

# We use Cramer's V first, since it takes less clusters
signif_df <- subset(results_df, cramersV > 0.15 & FDR < 0.01)
strong_signif_df <- subset(results_df, cramersV > 0.20 & FDR < 0.01)

# Plot again
par(mfrow = c(1,2))
p2 <- signif_df$FDR
p2[p2 == 0] <- min(signif_df$FDR[signif_df$FDR != 0]) # the true minimum
hist(-log10(p2),
     main = "FDR distribution\n(Cramer's V > 0.15)",
     xlab = "FDR",
     col = "skyblue",
     border = "white",
     breaks = 30)

p3 <- strong_signif_df$FDR
p3[p3 == 0] <- min(strong_signif_df$FDR[strong_signif_df$FDR != 0]) # the true minimum
hist(-log10(p3),
     main = "FDR distribution\n(Cramer's V > 0.20)",
     xlab = "FDR",
     col = "skyblue",
     border = "white",
     breaks = 30)
par(mfrow = c(1,1))

# Clusters closer in statistics to pipolb
sum(
  signif_df$presence >= (sum(clusters$pipolb) * 0.5) 
  & signif_df$presence <= (sum(clusters$pipolb) * 1.5)
)

sum(
  signif_df$percentage >= ((sum(clusters$pipolb)/nrow(clusters)*100) * 0.5) 
  & signif_df$percentage <= ((sum(clusters$pipolb)/nrow(clusters)*100) * 1.5)
)

#=========================================================================
# Step 4 - Extract data
#=========================================================================
#write.table(reults_df, file = "panta_clusters/ecoli_all_clusters.txt", sep = "\t", quote = F, row.names = F)

#write.table(signif_df$gene, file = "panta_clusters/ecoli_significant_clusters.txt", sep = "\t", quote = F, row.names = F, col.names = F)
#write.table(signif_df, file = "panta_clusters/ecoli_significant_clusters_dataset.txt", sep = "\t", quote = F, row.names = F,)

#writeLines(ecoli_dataset_phylogroup$assembly[ecoli_dataset_phylogroup$pipolin_presence == T], con = "../scripts_python/ecoli_pipolin_genomes.txt")
