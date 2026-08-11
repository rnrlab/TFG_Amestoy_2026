#### Analysis of relationship between clusters and phylogeny ####

"
We work with the results of functional annotation of clusters
with HHblits, combined with the phylogeny extracted from FastTree.

The goal is to see if the phylogeny relates in an way with the co-ocurrence
of certain clusters and the piPolB.
"

library(ape)
library(ggplot2)
library(parallel)
library(patchwork)

BiocManager::install("ggtree")
library(ggtree)
library(treeio)
library(dplyr)


#=========================================================================
# Step 0 - INPUTS
#=========================================================================
# Tree
tree <- read.tree("panta_clusters/tree")

# Matrix
core_clusters <- read.table(file = "panta_clusters/t_20_core_gene_presence_absence.Rtab", sep = "\t", header = T)

filtered_clusters <- read.table(file = "panta_clusters/out_pipolin_matrix.txt", sep = "\t", header = T)

# Data
ecoli_dataset_phylogroup <- read.table(file = 'ecoli_dataset_phylogroups.txt', sep = '\t', header = T)

signif_out_pipolin <- read.table(file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", header = T)


# Create piPolB presence-absence vector
filtered_clusters$pipolb <- ecoli_dataset_phylogroup$pipolin_presence[match(filtered_clusters$assembly, ecoli_dataset_phylogroup$assembly)]
sum(filtered_clusters$pipolb == 1) == sum(ecoli_dataset_phylogroup$pipolin_presence == 1)

#=========================================================================
# Step 1 - VIUALIZE TREE
#=========================================================================
# Check
head(tree$tip.label)
head(filtered_clusters$assembly)
sum(tree$tip.label %in% filtered_clusters$assembly)

# Visualize piPolB
p <- ggtree(tree, layout = "circular") %<+% filtered_clusters +
  geom_tippoint(aes(color = as.factor(pipolb)), size = 0.8) +
  scale_color_manual(values = c("0" = "#56B4E9", "1" = "#E41A1C")) +
  theme_tree()
ggsave("figuras/arbol_circular_pipolb.png", plot = p, width = 20, height = 20)

# Visualize MADS with less pYcY
p2 <- ggtree(tree, layout = "circular") %<+% filtered_clusters +
  geom_tippoint(aes(color = as.factor(groups_8688)), size = 0.8) +
  scale_color_manual(values = c("0" = "#56B4E9", "1" = "#EE5C42")) +
  theme_tree()
ggsave("figuras/arbol_circular_groups_8688.pdf", plot = p2, width = 20, height = 20)

pcomb <- p + p2
ggsave("figuras/arbol_circular_pipolb_groups_8688.pdf", plot = pcomb, width = 20, height = 20)

# Visualize cluster with 1.01% Prevalence
p3 <- ggtree(tree, layout = "circular") %<+% filtered_clusters +
  geom_tippoint(aes(color = as.factor(groups_25796)), size = 0.8) +
  scale_color_manual(values = c("0" = "#56B4E9", "1" = "green")) +
  theme_tree()
ggsave("figuras/arbol_circular_groups_25796.pdf", plot = p3, width = 20, height = 20)

#=========================================================================
# Step 2 - Statistical analysis
#=========================================================================
### A) Prepare data 
# Significative clusters outside pipolin
genes <- colnames(filtered_clusters)[-1]

# 20 clusters from the core genome
colnames(core_clusters)[colnames(core_clusters) == "X"] <- "assembly"
core_genes <- colnames(core_clusters)[-1]


tips_all <- tree$tip.label

### B) Moran's I + k-NN functions
# Calculation
calc_Moran_knn_fast <- function(tree, df, gene, tips, k = 5) {
  
  sub_tree <- drop.tip(tree, setdiff(tree$tip.label, tips))
  D <- cophenetic(sub_tree)
  
  sub_df <- df[df$assembly %in% tips, ]
  sub_df <- sub_df[match(tips, sub_df$assembly), ]
  
  x <- sub_df[[gene]]
  names(x) <- tips
  
  x <- x - mean(x, na.rm = TRUE)
  
  n <- length(tips)
  
  nn_idx <- t(apply(D, 1, order))[ , 2:(k+1), drop = FALSE]
  
  Wx <- numeric(n)
  
  for (i in 1:n) {
    Wx[i] <- sum(x[nn_idx[i, ]])
  }
  
  num <- sum(x * Wx)
  denom <- sum(x^2)
  
  (n / (n * k)) * num / denom
}

# Wrapper
calc_I_obs <- function(tree, df, gene,
                       reps = 50,
                       n_subsample = 300,
                       k = 5) {
  
  Is <- replicate(reps, {
    
    tips <- sample(tree$tip.label, n_subsample)
    
    tryCatch(
      calc_Moran_knn_fast(tree, df, gene, tips, k),
      error = function(e) NA
    )
  })
  
  Is <- Is[!is.na(Is)]
  
  mean(Is)
}

# Permutation test
calc_Moran_perm <- function(tree, df, gene,
                            reps = 500,
                            n_subsample = 300,
                            k = 5) {
  
  # observed
  I_obs <- calc_I_obs(tree, df, gene,
                      reps = 50,
                      n_subsample = n_subsample,
                      k = k)
  
  # null distribution
  I_null <- replicate(reps, {
    
    df_perm <- df
    df_perm[[gene]] <- sample(df_perm[[gene]])
    
    calc_I_obs(tree, df_perm, gene,
               reps = 10,
               n_subsample = n_subsample,
               k = k)
  })
  
  # p-value por permutacion --> mínimo pvalue = 1/(n+1)
  p_value <- (sum(abs(I_null) >= abs(I_obs)) + 1) /
    (length(I_null) + 1)
  
  list(
    I_obs = I_obs,
    p_value = p_value,
    null_mean = mean(I_null),
    null_sd = sd(I_null)
  )
}
calc_Moran_perm(tree = tree, df = core_clusters, gene = "groups_621", reps = 50, n_subsample = 200, k = 5)


### B) Moran's I + k-NN main loop
morans_out <- mclapply(core_genes, function(gene){
  
  calc_Moran_perm(
    tree,
    core_clusters,
    gene,
    reps = 500,
    n_subsample = 300,
    k = 5
  )
  
}, mc.cores = detectCores() - 1)

#=========================================================================
# Step 3 - Parse results to dataframe
#=========================================================================

morans_df <- data.frame(
  gene = core_genes,
  I_obs = sapply(morans_out, `[[`, "I_obs"),
  p_value = sapply(morans_out, `[[`, "p_value"),
  null_mean = sapply(morans_out, `[[`, "null_mean"),
  null_sd = sapply(morans_out, `[[`, "null_sd")
)

# FDR correction
morans_df$FDR <- p.adjust(morans_df$p_value, method = "BH")

# Classification
morans_df$phylo_status <- ifelse(
  morans_df$FDR < 0.05,
  "non_random",
  "random"
)

# Output table
write.table(
  morans_df,
  "panta_clusters/20_core_morans_knn_perm_results.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#####################
#m <- read.table(file = "panta_clusters/20_core_morans_knn_perm_results.txt", sep = "\t", header = T)
