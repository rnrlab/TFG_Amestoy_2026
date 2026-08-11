#### Analysis of IntegronFinder2.0 results over 16417 genomes with pipolins ####

"
We use the results from integron_finder (IntegronFinder2.0) metadata as input.

We work with 2 main files:
 - pipolin_metadata_filtered2.txt --> summary table of all MGEs found on bacterial genomes with pipolins
 - info_from_integronfinder.txt --> summary table of all Integrons found on bacterial genomes with pipolins
"

packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "readr")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=========================================================================
#===== Step 0: Importing data ============================================
#=========================================================================

## A) Import
# Genomes w/ pipolins
pipolin_metadata_filtered2 <- read.table(
  file = "pipolin_metadata_filtered2.txt",
  sep = '\t',
  header = TRUE
)

# Genomes w/ pipolins and integrons
integrons_per_genome <- read.table(
  file = "info_from_integronfinder.txt",
  sep = '\t',
  header = TRUE
)


## B) Merge and process
# Left-join --> only keep genomes with pipolins
pipolin_metadata_filtered2 <- merge(pipolin_metadata_filtered2, integrons_per_genome, by = 'assembly', all.x = TRUE)

# Change "NA" for "0"
pipolin_metadata_filtered2[c('total_integrons', 'calin', 'complete_integrons', 'In0')][
  is.na(pipolin_metadata_filtered2[c('total_integrons', 'calin', 'complete_integrons', 'In0')])
] <- 0

#=========================================================================
#===== Step 1: Distribution of Integrons per taxon =======================
#=========================================================================

### A) Parse data
# Make long table
integrons_genus <- pipolin_metadata_filtered2 %>%
  group_by(class, order, genus) %>%
  summarise(
    total_genomes = n(),
    genomes_with_integrons = sum(complete_integrons > 0, na.rm = TRUE),
    complete_integrons = sum(complete_integrons, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_with_integrons = genomes_with_integrons / total_genomes * 100
  )
integrons_genus[c('class', 'order', 'genus')] <- lapply(integrons_genus[c('class', 'order', 'genus')], as.factor)

# Filter
integrons_genus <- integrons_genus %>%
  filter(complete_integrons > 0)

### B) Plot percentages by genus
ggplot(
  data = integrons_genus %>% 
    filter(genomes_with_integrons > 5), 
  aes(
    x = reorder(genus, pct_with_integrons), 
    y = pct_with_integrons, 
    fill = genus
  )
) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
  
  coord_flip() +
  
  scale_fill_manual(values = c(
    "Aeromonas" = "#FA8072",
    "Staphylococcus" = "#ffa947",
    "Citrobacter" = "#5bad51",
    "Enterobacter" = "#7FE072",
    "Escherichia" = "#25d9c1",
    "Limosilactobacillus" = "#4781e6",
    "Pseudosulfitobacter" = "#c782ff",
    "Vibrio" = "#e32bac",
    "Methylobacterium" = "#8a2b32",
    "Salmonella" = "#ff94f1",
    "Klebsiella" = "#1f7a8c"
  )) +
  
  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, 5)
  ) +
  
  theme_classic() +

  labs(
    title = "Porcentaje de genomas con integrones completos por género",
    x = "Género",
    y = "Integrones (%)",
    fill = "Género"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    legend.text = element_text(face = "italic", size = 18),
    strip.text = element_text(size = 22, face = "italic"),
    plot.title = element_text(size = 25, face = "bold", vjust = 2),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23, face = "italic")
  )


#=========================================================================
#===== Step 2: Distribution of Integrons in Escherichia coli =============
#=========================================================================



#=========================================================================
#===== Step 3: Exporting data ============================================
#=========================================================================
# Save big dataset to file
write.table(
  pipolin_metadata_filtered2,
  file = 'pipolin_metadata_filtered3.txt',
  sep = '\t',
  row.names = FALSE
)
