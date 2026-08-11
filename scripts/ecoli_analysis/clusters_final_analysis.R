#### Final analysis of clusters obtained by PanTA ####

"
Here we import the functional-annotation results obtained both
by PADLOC and DeepKOALA (for defense systems and general annotation),
and combine them with the metadata of the Pangenome clusters
with a significant correlation with the piPolB gene,
given that the proteins within these clusers are outside the pipolin
(by 'given' I mean we made sure of it).
"

packages <- c("dplyr", "ggplot2", "stringr", "patchwork")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
# Step 0 - INPUTS
#=========================================================================
all_clusters <- read.table(file = "panta_clusters/ecoli_all_clusters.txt", sep = "\t", header = T)

signif_out_pipolin <- read.table(file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", header = T)

padloc <- read.csv(file = "panta_clusters/padloc_best_hits_per_cluster.csv", header = T)

# All DeepKOALA results
koala_all <- read.csv("panta_clusters/deepkoala_all_results.csv", header = TRUE)

#=========================================================================
# Step 1 - Parse PADLOC results
#=========================================================================
### A) Add to complete dataframe
padloc <- padloc %>%
  rename(
    gene = cluster,
    padloc_best_hmm = best_hmm,
    padloc_best_target = best_target,
    padloc_best_evalue = best_evalue,
    padloc_best_score = best_score,
    padloc_bias = bias,
    padloc_hmm_name = hmm_name,
    padloc_hmm_accession = hmm_accession,
    padloc_hmm_description = hmm_description,
    padloc_hmm_length = hmm_length,
    padloc_system = system
  )

#signif_out_pipolin <- merge(signif_out_pipolin, padloc, by = "gene", all.x = TRUE)
#signif_out_pipolin$defense_system <- ifelse(!is.na(signif_out_pipolin$padloc_best_hmm), 1, 0)
#signif_out_pipolin$correlation <- ifelse(signif_out_pipolin$phi >= 0, "positive", "negative")

### B) Filter by defense system presence
signif_out_pipolin_DS <- signif_out_pipolin[signif_out_pipolin$defense_system == T,]

# Filter by e-value < 0.01
summary(signif_out_pipolin_DS$padloc_best_evalue)
signif_out_pipolin_DS <- signif_out_pipolin_DS[signif_out_pipolin_DS$padloc_best_evalue < 1e-03,]

### C) Plot abundance of defense systems
sort(table(signif_out_pipolin_DS$padloc_system), decreasing = TRUE)

# Frequency dataframe
DS_freq_df <- signif_out_pipolin_DS %>%
  count(padloc_system, correlation) %>% 
  group_by(padloc_system) %>%
  mutate(
    total = sum(n),
    pos = ifelse(correlation == "positive", n, 0),
    neg = ifelse(correlation == "negative", n, 0),
    score = (sum(pos) - sum(neg)) / total,  
    pct = 100 * n / total
  ) %>%
  ungroup()


# Correct some groups
DS_freq_df <- bind_rows(
  DS_freq_df,
  data.frame(padloc_system = "Cas", 
             correlation = "positive", 
             n = 3,
             total = 3,
             pos = 3,
             neg = 0,
             score = 1.000,
             pct = 100.00)
)
DS_freq_df <- DS_freq_df[!(DS_freq_df$padloc_system %in% c("casc", "casd", "case")),]

# Plot
ds_plot <- ggplot(DS_freq_df, aes(x = n, y = reorder(padloc_system, total), fill = correlation)) +
  geom_col(position = "stack") +
  theme_bw() +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  geom_text(
    aes(x = total, label = paste0("S: ", round(score, 2))),
    hjust = -0.1,
    size = 4,
    fontface = "bold"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    title = "Frecuencia de genes de defensa",
    x = "Nº clusters",
    y = "Gen",
    fill = "Correlación"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
ds_plot

#=========================================================================
# Step 2 - Parse DeepKOALA + Kegg Orthology results
#=========================================================================
### A) Filter: Above threshold --> if multiple | if not --> ↑Probability
best_koala <- koala_all %>%
  group_by(name) %>%
  arrange(desc(annotate == "*"), desc(probability)) %>% 
  slice(1) %>% 
  ungroup()
sum(best_koala$annotate == "*")

### B) Kegg Mapper IDs
length(unique(best_koala$predict_label))
(kegg_ids <- unique(best_koala$predict_label))

# Export Kegg IDs as INPUT --> use Kegg Mapper to extract function
writeLines(kegg_ids, con = "panta_clusters/input_ids_for_kegg_mapper.txt")

# Import Kegg Mapper OUTPUT
kegg_annotations <- read.csv2("panta_clusters/ko.txt", header = T)

# Check if KO classes are the same
all(best_koala$predict_label %in% kegg_annotations$reaction_ko)

# Merge Koala dataset with KO classes info
best_koala <- merge(best_koala, kegg_annotations, by.x = "predict_label", by.y = "reaction_ko")
length(unique(best_koala$name))
length(best_koala$predict_label)
length(unique(best_koala$predict_label))

best_koala <- best_koala %>%
  rename(gene = name)

### C) Plot KO categories
## C.1) Parse beforehand
# Add correlation variable
best_koala <- merge(best_koala, 
                    signif_out_pipolin[, c("gene", "correlation")], 
                    by = "gene")

# Add defense system presence variable
best_koala <- merge(best_koala, 
                    signif_out_pipolin[, c("gene", "defense_system")], 
                    by = "gene")

# Add combined correlation X defense system variable
best_koala$correlationXdefense <- with(best_koala,
  ifelse(correlation == "positive" & defense_system == TRUE,  "C+D+",
  ifelse(correlation == "positive" & defense_system == FALSE, "C+D-",
  ifelse(correlation == "negative" & defense_system == TRUE,  "C-D+", "C-D-")))
)

# Add "mads" presence
signif_out_pipolin$mads_presence <- ifelse(
  !is.na(signif_out_pipolin$padloc_system) & signif_out_pipolin$padloc_system == "mads",
  "mads+",
  "mads-"
)

best_koala <- merge(best_koala, 
                    signif_out_pipolin[, c("gene", "mads_presence")], 
                    by = "gene")

# Add combined correlation X mads_presence variable
best_koala$correlationXmads <- with(best_koala,
  ifelse(correlation == "positive" & mads_presence == "mads+",  "C+mads+",
  ifelse(correlation == "positive" & mads_presence == "mads-", "C+mads-",
  ifelse(correlation == "negative" & mads_presence == "mads+",  "C-mads+", "C-mads-")))
)


## C.2) Create counts xtable
# Count correlation per pathways per class combination
xkegg <- best_koala %>%
  count(a_class, b_class, pathway, correlationXdefense, name = "Freq")
table(xkegg$Freq)

## C.3) Plot
# Correlation X Defense system
ggplot(xkegg %>% filter(Freq >= 5), aes(y = Freq, x = pathway, fill = correlationXdefense)) + 
  geom_col(
    position = "stack",
    color = "grey40",
    linewidth = 0.2,
    alpha = 0.8
  ) +
  facet_wrap(~ a_class + b_class, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "C+D+" = "#A52A2A",
    "C+D-" = "#EE5C42",
    "C-D+" = "#2166ac",
    "C-D-" = "#67a9cf"
  ))

# Correlation X Mads
xkegg2 <- best_koala %>%
  count(a_class, b_class, pathway, correlationXmads, name = "Freq")
table(xkegg$Freq)

ggplot(xkegg2 %>% filter(Freq >= 7), aes(y = Freq, x = pathway, fill = correlationXmads)) + 
  geom_col(
    position = "stack",
    color = "grey40",
    linewidth = 0.2,
    alpha = 0.8
  ) +
  facet_wrap(~ a_class + b_class, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "C+mads+" = "#A52A2A",
    "C+mads-" = "#EE5C42",
    "C-mads+" = "#2166ac",
    "C-mads-" = "#67a9cf"
  ))

# Both Correlation X Defense system X mads presence w/ facet_wrap
# We know that "mads" is only present in C-D+
xkegg3 <- best_koala %>%
  count(a_class, b_class, pathway, correlationXdefense, mads_presence, name = "Freq") %>%
  mutate(fill_group = ifelse(correlationXdefense == "C-D+" & mads_presence == "mads+",
                             "C-D+ mads+",
                             correlationXdefense))
KO_corr_def_mads_1 <- ggplot(
  data = xkegg3 %>% 
    filter(Freq >= 5) %>% 
    filter(b_class %in% c("09101 Carbohydrate metabolism", "09131 Membrane transport", "09132 Signal transduction")),
  aes(y = Freq, x = pathway, fill = fill_group)) + 
  geom_col(
    position = "stack",
    color = "grey40",
    linewidth = 0.2,
    alpha = 0.8
  ) +
  facet_grid(~ a_class + b_class, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c(
    "C+D+" = "#A52A2A",
    "C+D-" = "#EE5C42",
    "C-D+" = "#2166ac",
    "C-D-" = "#67a9cf",
    "C-D+ mads+" = "#08306B"
  )) +
  labs(
    title = "Metabolic categorization",
    x = "",
    y = "Frecuencia",
    fill = "Correlación-Defensa"
  ) +
  theme(
    axis.title.x = element_text(size = 25),
    #axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    title = element_text(size = 25),
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 17, face = "plain")
  ) +
  scale_y_continuous(limits = c(0,60))
KO_corr_def_mads_1

KO_corr_def_mads_2 <- ggplot(
  data = xkegg3 %>% 
    filter(Freq >= 5) %>% 
    filter(b_class %in% c("09182 Protein families: genetic information processing", "09183 Protein families: signaling and cellular processes")),
  aes(y = Freq, x = pathway, fill = fill_group)) + 
  geom_col(
    position = "stack",
    color = "grey40",
    linewidth = 0.2,
    alpha = 0.8
  ) +
  facet_grid(~ a_class + b_class, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c(
    "C+D+" = "#A52A2A",
    "C+D-" = "#EE5C42",
    "C-D+" = "#2166ac",
    "C-D-" = "#67a9cf",
    "C-D+ mads+" = "#08306B"
  )) +
  labs(
    title = "",
    x = "",
    y = "Frecuencia",
    fill = "Correlation-Defense"
  ) +
  theme(
    axis.title.x = element_text(size = 25),
    #axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    title = element_text(size = 25),
  ) +
  theme(
    legend.position = "",
    legend.title = element_text(size = 10),
    strip.text = element_text(size = 21, face = "plain")
  ) +
  scale_y_continuous(limits = c(0,60))
KO_corr_def_mads_2

#=========================================================================
# Step 3 - FDR, Phi and Moran's I distribution
#=========================================================================
### A) All clusters
# Chi-Square
ggplot(all_clusters, aes(x = chi.s)) +
  geom_histogram(position = "identity", bins = 90, fill = "#67a9cf", color = "black") +
  theme_minimal() +
  labs(title = "Distribución de valores de Chi-cuadrado", x = "Chi-Cuadrado", y = "Frecuencia")

# FDR - correlation
pseudo <- min(all_clusters$FDR.corr[all_clusters$FDR.corr > 0]) / 10

all_clusters$logFDR <- -log10(
  ifelse(
    all_clusters$FDR.corr == 0,
    pseudo,
    all_clusters$FDR.corr)
  )
ggplot(all_clusters, aes(x = logFDR)) +
  geom_histogram(position = "identity", bins = 50, fill = "#67a9cf", color = "black") +
  geom_vline(
    xintercept = -log(0.001),
    linetype = "dashed",
    linewidth = 2,
    color = "red"
  ) +
  annotate("text", x = -log(0.001), y = Inf,
           label = "FDR = 0,001",
           vjust = 2, hjust = -0.05, size = 15) +
  labs(title = "Distribución de FDR (Chi-2)", 
       x = expression(-log[10](FDR)), 
       y = "Frecuencia") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(size = 33),
    axis.text.y = element_text(size = 33),
    plot.title = element_text(size = 30, face = "bold", vjust = 2)
  )

# Cramer's V
ggplot(all_clusters, aes(x = cramersV)) +
  geom_histogram(position = "identity", bins = 75, fill = "#67a9cf", color = "black") +
  geom_vline(
    xintercept = 0.15,
    linetype = "dashed",
    linewidth = 2,
    color = "red"
  ) +
  annotate("text", x = 0.15, y = Inf,
           label = "V = 0,15",
           vjust = 2, hjust = -0.05, size = 15) +
  labs(title = "Distribución de V de Cramer", 
       x = "V de Cramer", 
       y = "Frecuencia") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 27, face = "bold"),
    axis.title.y = element_text(size = 27, face = "bold"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    plot.title = element_text(size = 26, face = "bold", vjust = 2)
  )

### B) Significant clusters outside the pipolin
## Cramer's V
# Boxplot
ggplot(signif_out_pipolin, aes(x = correlation, y = cramersV, fill = correlation)) +
  geom_boxplot() +
  scale_fill_manual(values = c("negative" = "#8DA0CB", "positive" = "#66C2A5")) +
  theme_minimal() +
  labs(
    title = "Distribución de Cramer's V",
    x = NULL,
    y = "Cramer's V") +
  coord_flip() +
  scale_y_continuous(
    limits = c(0.14, 0.40), 
    breaks = seq(0, 1, by = 0.03)) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

# Histogram
(pos_median <- median(signif_out_pipolin$cramersV[signif_out_pipolin$correlation == "positive"]))
(neg_median <- median(signif_out_pipolin$cramersV[signif_out_pipolin$correlation == "negative"]))
median_df <- data.frame(
  correlation = c("positive", "negative"),
  median = c(pos_median, neg_median)
)

ggplot(signif_out_pipolin, aes(x = cramersV)) +
  geom_histogram(position = "identity", bins = 20, fill = "#67a9cf", color = "black") +
  geom_vline(
    data = median_df,
    aes(xintercept = median),
    linetype = "dashed",
    linewidth = 1,
    color = "black"
  ) +
  geom_text(
    data = median_df,
    aes(
      x = median,
      y = Inf,
      label = paste0("Mediana = ", round(median, 2))
    ),
    vjust = 2,
    hjust = -0.11,
    size = 6
  ) +
  facet_wrap(. ~ correlation) +
  labs(title = "Distribución de valores de la V de Cramer", 
       x = "V de Cramer", 
       y = "Frecuencia") +
  scale_x_continuous(breaks = seq(0.15, 0.5, 0.05)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 20),
    strip.text = element_text(size = 18, face = "bold")
  )

# Dispersion: -log(FDR) x Phi 
pseudo_signif_out <- min(signif_out_pipolin$FDR.corr[signif_out_pipolin$FDR.corr > 0]) / 10

signif_out_pipolin$logFDR <- -log10(
  ifelse(
    signif_out_pipolin$FDR.corr == 0,
    pseudo_signif_out,
    signif_out_pipolin$FDR.corr)
)

signif_out_pipolin <- signif_out_pipolin %>%
  relocate(logFDR, .after = FDR.corr)

ggplot(signif_out_pipolin, aes(x = logFDR, y = phi, color = correlation)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_hline(
    data = median_df,
    aes(yintercept = 0),
    linetype = "dashed",
    linewidth = 1.5,
    color = "black"
  ) +
  scale_color_manual(
    breaks = c("positive", "negative"),
    values = c(
      "positive" = "#A52A2A",
      "negative" = "#67a9cf"),
    labels = c(
      "positive" = "Positiva",
      "negative" = "Negativa"),
    name = "Correlación"
  ) +
  labs(title = "Distribución de Phi según correlación", 
       x = expression(-log[10](FDR)), 
       y = "Phi") +
  theme_classic(base_size = 16) +

  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = 2),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 33),
    axis.text.y = element_text(size = 33),
    legend.text = element_text(size = 33),
    legend.title = element_text(size = 33)
  )

nrow(signif_out_pipolin[signif_out_pipolin$correlation == "positive",])
nrow(signif_out_pipolin[signif_out_pipolin$correlation == "negative",])





ggplot(signif_out_pipolin, aes(x = pYcY, y = phi, color = correlation)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "positive" = "#A52A2A",
      "negative" = "#67a9cf"),
    labels = c(
      "positive" = "Positiva",
      "negative" = "Negativa"),
    name = "Correlación"
  ) +
  labs(title = "Distribución de valores de Phi en función de su significancia estadística", 
       x = expression(-log[10](FDR)), 
       y = "Phi") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 26)
  )


# Moran's I
ggplot(signif_out_pipolin, aes(x = correlation, y = morans_I_obs, fill = correlation)) +
  geom_boxplot() +
  scale_fill_manual(values = c("negative" = "#8DA0CB", "positive" = "#66C2A5")) +
  theme_minimal() +
  labs(
    title = "Distribución de Moran's I",
    x = NULL,
    y = "Moran's I") +
  coord_flip() +
  scale_y_continuous(
    limits = c(-0.05, 0.05)) +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

ggplot(signif_out_pipolin, aes(x = morans_I_obs)) +
  geom_histogram(position = "identity", bins = 40, fill = "#67a9cf", color = "black") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 2,
    color = "black") +
  labs(title = "Distribución de valores del Índice de Moran", 
       x = "Índice de Moran", 
       y = "Frecuencia") +
  scale_x_continuous(
    limits = c(-0.03, 0.03),
    breaks = seq(-0.03, 0.03, 0.01)) +
  theme_classic() +
  
  theme(
    plot.title = element_text(size = 30, face = "bold", vjust = 2),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 33),
    axis.text.y = element_text(size = 33),
    legend.text = element_text(size = 33),
    legend.title = element_text(size = 33)
  )

# P.value - Moran's I
ggplot(signif_out_pipolin, aes(x = correlation, y = morans_p_value, fill = correlation)) +
  geom_boxplot() +
  scale_fill_manual(values = c("negative" = "#8DA0CB", "positive" = "#66C2A5")) +
  theme_minimal() +
  labs(
    title = "Distribución de Moran's p-value",
    x = NULL,
    y = "p-value") +
  coord_flip() +
  theme(
    legend.position = 'right',
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

#=========================================================================
# Step 4 - Genes annotation frecuency
#=========================================================================
### A) Create more general functional annotation and count
#signif_out_pipolin <- signif_out_pipolin %>% mutate(name_hh_group = str_split(name_hh, "[-_]") %>% sapply(`[`, 1)) %>% relocate(name_hh_group, .after = name_hh)
#signif_out_pipolin_DS <- signif_out_pipolin_DS %>% mutate(name_hh_group = str_split(name_hh, "[-_]") %>% sapply(`[`, 1)) %>% relocate(name_hh_group, .after = name_hh)

out_clusters_counts <-  signif_out_pipolin %>%
  count(name_hh_group, correlation)

### B) Relación HHblits - PADLOC
out_clusters_counts$name_hh_group[out_clusters_counts$name_hh_group %in% signif_out_pipolin_DS$name_hh_group]

padloc_map <- signif_out_pipolin_DS %>%
  group_by(name_hh_group) %>%
  summarise(
    padloc_system = paste(unique(padloc_system), collapse = "; "),
    .groups = "drop"
  )

out_clusters_counts <- out_clusters_counts %>%
  left_join(padloc_map, by = "name_hh_group")

out_clusters_counts$final_name <- ifelse(
  is.na(out_clusters_counts$padloc_system),
  out_clusters_counts$name_hh_group,
  out_clusters_counts$padloc_system
)

# Rename some clusters
out_clusters_counts$final_name[out_clusters_counts$final_name == "retron; pd-t7-2"] <- "pd-t7-2"
out_clusters_counts$final_name[out_clusters_counts$final_name == "casc; casd"] <- "Cas"
out_clusters_counts$final_name[out_clusters_counts$final_name == "case"] <- "Cas"

# New table with corrected groups
corrected_clusters_counts <- out_clusters_counts %>%
  group_by(final_name, correlation) %>%
  summarise(
    total_n = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(correlation, desc(total_n))

# Correct some clusters counts
corrected_clusters_counts$total_n[
  corrected_clusters_counts$final_name == "PTS"
] <- 5
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("PTSIIB")),]

corrected_clusters_counts$total_n[
  corrected_clusters_counts$final_name == "Fimbrial"
] <- 8

# Add some other interesting clusters
corrected_clusters_counts <- bind_rows(
  corrected_clusters_counts,
  data.frame(final_name = "IS-Transposase", correlation = "positive", total_n = 3)
)
corrected_clusters_counts <- bind_rows(
  corrected_clusters_counts,
  data.frame(final_name = "IS-Transposase", correlation = "negative", total_n = 2)
)
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("DDE")),]

corrected_clusters_counts <- bind_rows(
  corrected_clusters_counts,
  data.frame(final_name = "Tox-Antitox", correlation = "positive", total_n = 10)
)
corrected_clusters_counts <- bind_rows(
  corrected_clusters_counts,
  data.frame(final_name = "Tox-Antitox", correlation = "negative", total_n = 4)
)
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("Toxin")),]
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("Tox")),]
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("Conotoxin")),]
corrected_clusters_counts <- corrected_clusters_counts[!(corrected_clusters_counts$final_name %in% c("Ntox34")),]

# Stablisth defense system status
defense_names <- out_clusters_counts %>%
  filter(!is.na(padloc_system)) %>%
  pull(final_name) %>%
  unique()

corrected_clusters_counts <- corrected_clusters_counts %>%
  mutate(
    defense_system = final_name %in% defense_names
  )

# Filter small groups by Cramer's V
lookup_cramers <- signif_out_pipolin %>%
  mutate(
    lookup_name = coalesce(padloc_system, name_hh_group)
  ) %>%
  group_by(lookup_name) %>%
  summarise(
    cramersV = max(cramersV),
    .groups = "drop"
  )

corrected_clusters_counts_f <- corrected_clusters_counts %>%
  left_join(
    lookup_cramers,
    by = c("final_name" = "lookup_name")
  ) %>%
  mutate(
    cramersV_small = ifelse(total_n == 1, cramersV, NA)
  ) %>%
  select(-cramersV)

corrected_clusters_counts_f <- corrected_clusters_counts_f %>%
  filter(
    total_n > 1 |
      (total_n == 1 & !is.na(cramersV_small))
  )

corrected_clusters_counts_f <- corrected_clusters_counts_f %>%
  filter(
    total_n > 1 |
      (total_n == 1 & !is.na(cramersV_small) & cramersV_small >= 0.15)
  )

# Rename some clusters
corrected_clusters_counts_f <- corrected_clusters_counts_f[!(corrected_clusters_counts_f$final_name %in% c("F")),]

corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Peripla"] <- "Peripla-BP"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Ni"] <- "Ni-hydr-CYTB"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "MtrB"] <- "MtrB-PioB"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "BPD"] <- "BPD-transp-2"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Aldo"] <- "Aldo-ket-red"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Zn"] <- "Zn-Tnp-IS1"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "PapC"] <- "PapC-N"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "LysR"] <- "LysR-substrate"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "LT"] <- "LT-IIB"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Glu"] <- "Glu-dehyd-C"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "EF"] <- "EF-hand-3"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "CSN4"] <- "CSN4-RPN5-eIF3a"
corrected_clusters_counts_f$final_name[corrected_clusters_counts_f$final_name == "Abhydrolase"] <- "Abhydrolase-3"
corrected_clusters_counts_f <- corrected_clusters_counts_f[!(corrected_clusters_counts_f$final_name %in% c("MqsR")),]

### C) Plot frequencies
table(corrected_clusters_counts_f$total_n)
ggplot(corrected_clusters_counts_f %>%
         filter(total_n > 1 | (total_n == 1 & !is.na(cramersV_small) & cramersV_small > 0.1925203)), 
       aes(x = reorder(final_name, total_n), y = total_n, fill = correlation)) +
  geom_col(position = "dodge", color = "black", width = 0.7) +
  labs(
       title = "Frecuencia de genes que correlacionan con la piPolB", 
       x = "", 
       y = "Nº Clusters",
       fill = "Correlación") +
  scale_fill_manual(
    values = c(
      "negative" = "#67a9cf",
      "positive" = "#A52A2A"),
    labels = c(
      "negative" = "Negativa",
      "positive" = "Positiva"
    )) +
  scale_y_continuous(breaks = seq(0, 30, 2)) +
  facet_wrap(. ~ defense_system, scales = "free") +
  facet_wrap(~ defense_system, scales = "free",
             labeller = labeller(defense_system = c(
               "TRUE" = "Genes de defensa",
               "FALSE" = "Anotación funcional"
             ))) +
  coord_flip() +
  theme_classic() +

  theme(
    plot.title = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 21),
    strip.text = element_text(size = 18, face = "bold")
  )

### C) Plot frequencies filtered (for better visualization)
ggplot(corrected_clusters_counts_f %>%
         filter(final_name %in% c(
           "Fimbrial", "Tox-Antitox", "T6SS", "PTS",
           "SIS", "Peripla-BP", "IS-Transposase",
           "mads", "pd-t7-2", "Cas", "lamassu-fam", "dsr1"
         )), 
       aes(x = reorder(final_name, total_n), y = total_n, fill = correlation)) +
  geom_col(position = "dodge", color = "black", width = 0.7) +
  labs(
    title = "Frecuencia de genes que correlacionan con la piPolB", 
    x = "", 
    y = "Nº Clústeres",
    fill = "Correlación") +
  scale_fill_manual(
    values = c(
      "negative" = "#67a9cf",
      "positive" = "#A52A2A"),
    labels = c(
      "negative" = "Negativa",
      "positive" = "Positiva"
    )) +
  scale_y_continuous(
    breaks = seq(0, 30, 4),
  ) +
  facet_wrap(. ~ defense_system, scales = "free",
             labeller = labeller(defense_system = c(
               "TRUE" = "Genes de defensa",
               "FALSE" = "Anotación funcional"
             ))) +
  coord_flip() +
  theme_classic() +
  
  theme(
    plot.title = element_text(size = 35, face = "bold"),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    axis.text.x = element_text(size = 33),
    axis.text.y = element_text(size = 45),
    legend.text = element_text(size = 39),
    legend.title = element_text(size = 39),
    legend.position = "",
    strip.text = element_text(size = 39, face = "bold"),
    plot.margin = margin(10, 35, 10, 10)
  )
ggsave(filename = "../escrito_TFG/figuras/genes.png", width = 26, height = 14)

#=========================================================================
# Step 5 - Export dataset with new columns
#=========================================================================
#write.table(signif_out_pipolin, file = "panta_clusters/ecoli_signif_out_pipolin_dataset.txt", sep = "\t", quote = F, row.names = F)

sum(signif_out_pipolin$phi < 0)
sum(signif_out_pipolin$phi > 0)
