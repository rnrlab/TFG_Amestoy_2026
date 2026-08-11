#### Analysis of relationship between clusters and phylogeny, Part 2 ####

"
We work with the results of the statistical analysis with the Moran's Index,
or Moran's I.
This statistics was obtained with the script 'tree_analysis.R'.
"
packages <- c("dplyr")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
# Step 0 - INPUTS
#=========================================================================
signif_out_pipolin <- read.table(file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", header = T)

morans_results <- read.table(file = "panta_clusters/morans_knn_perm_results.txt", sep = "\t", header = T)

#=========================================================================
# Step 1 - Parse results and compare
#=========================================================================
summary(morans_results)


# Merge
signif_out_pipolin <- merge(signif_out_pipolin, morans_results, by = "gene")

library(dplyr)

signif_out_pipolin <- signif_out_pipolin %>% 
  rename(
    FDR.corr = FDR.x,
    morans_I_obs = I_obs,
    morans_p_value = p_value,
    morans_null_mean = null_mean,
    morans_null_sd = null_sd,
    morans_FDR = FDR.y
  )

# Export
write.table(signif_out_pipolin, file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", quote = F, row.names = F)
