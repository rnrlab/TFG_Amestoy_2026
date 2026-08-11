#### Analysis of ClermonTyping and MLST results over 32949 E. coli genomes ####

"
We use as inputs the tabulted results of ClermonTyping (phylogroups)
and MLST (sequence types)

We work with files:
  - all_phylogroups.txt: tabulated text file with 1 row per genome with info about its chosen phylogroups.
  - pipolin_metadata_escherichia.txt: tabulated text file with plenty info about Escherichia genomes with pipolins. 
  - results_mlst.txt: tabulated text faile with 1 row per genome with info about its chosen ST.
"

packages <- c("ggplot2", "dplyr", "tidytext", "patchwork")
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
all_phylogroups <- read.table(
  file = 'all_phylogroups.txt',
  sep = '\t',
  header = TRUE
)

ecoli_info_from_EP <- read.csv(
  file = "ecoli_info_from_EP.txt",
  sep = "\t"
)

"
pipolin_metadata_escherichia <- read.table(
  file = 'pipolin_metadata_escherichia.txt',
  sep = '\t',
  header = TRUE
)
"

results_mlst <- read.table(
  file = 'results_mlst.txt',
  sep = '\t',
  header = TRUE
)


### B) Parse ClermonTyping data
## Groups in phylogroups
unique(all_phylogroups$phylogroup)

## Groups in mash_groups
unique(all_phylogroups$mash_group)

## Groups in phylogroups but not in mash_groups
setdiff(all_phylogroups$phylogroup, all_phylogroups$mash_group)

## Groups in mash_groups but not in phylogroups
setdiff(all_phylogroups$mash_group, all_phylogroups$phylogroup)

## Assemblies with difference between phylogroups and mash_groups
length(all_phylogroups$assembly[all_phylogroups$phylogroup != all_phylogroups$mash_group])
length(all_phylogroups$assembly[all_phylogroups$phylogroup == all_phylogroups$mash_group])
all_phylogroups$assembly[all_phylogroups$phylogroup != all_phylogroups$mash_group]

## Assemblies with only common phylogroups ("A","E","B2","B1","F","D","C","G")
# Only in phylogroups
length(all_phylogroups$assembly[all_phylogroups$phylogroup %in% c("A","E","B2","B1","F","D","C","G")])

# Only in mash_group
length(all_phylogroups$assembly[all_phylogroups$mash_group %in% c("A","E","B2","B1","F","D","C","G")])

# Phylogroups OR mash_group
length(all_phylogroups$assembly[all_phylogroups$phylogroup %in% c("A","E","B2","B1","F","D","C","G") | all_phylogroups$mash_group %in% c("A","E","B2","B1","F","D","C","G")])

# Phylogroups AND mash_group
length(all_phylogroups$assembly[all_phylogroups$phylogroup %in% c("A","E","B2","B1","F","D","C","G") & all_phylogroups$mash_group %in% c("A","E","B2","B1","F","D","C","G")])


### C) Parse ecoli_info_from_EP data
# Genomes with pipolins in our new ~33k E. coli genomes dataset
length(intersect(ecoli_info_from_EP$assembly, all_phylogroups$assembly)) # just checkin
nrow(ecoli_info_from_EP[ecoli_info_from_EP$pipolin_presence == 1,])


### D) Parse MLST data
length(unique(results_mlst$ST))
names(which.max(table(results_mlst$ST))) # more common ST


### E) Create combined dataset
ecoli_dataset <- ecoli_info_from_EP

# Add Phylogroups info
ecoli_dataset <- merge(ecoli_dataset, all_phylogroups, by = 'assembly')
length(ecoli_dataset$assembly[ecoli_dataset$pipolin_presence == 1]) # just checkin

# Add MLST info
ecoli_dataset <- merge(ecoli_dataset, results_mlst, by = 'assembly')
length(unique(ecoli_dataset$ST))

# Check percentage of genomes with pipolins
(nrow(ecoli_dataset[ecoli_dataset$pipolin_presence == 1,]) / nrow(ecoli_dataset)) * 100
(1198 / nrow(ecoli_dataset)) * 100

#=========================================================================
# Step 1: Plot phylogroups by pipolin presence
#=========================================================================
sapply(c("A", "B1", "B2", "C", "D", "E", "F", "G"), function(group){
  (sum(all_phylogroups[["phylogroup"]] == group) / 32700) * 100
})

# Make long table
ecoli_phylogroups_freq <- ecoli_dataset %>%
  group_by(phylogroup, pipolin_presence) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(phylogroup) %>%
  mutate(
    total = sum(count),
    pct_phylogroup = count / total * 100
  ) %>%
  ungroup()

# Plot
plot_phylogroups <- ggplot(ecoli_phylogroups_freq %>%
         filter(phylogroup %in% c("A","E","B2","B1","F","D","C","G")), 
       aes(x = phylogroup, y = count, fill = factor(pipolin_presence))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    data = ecoli_phylogroups_freq %>%
      filter(pipolin_presence == 1 & phylogroup %in% c("A","E","B2","B1","F","D","C","G")),
    aes(x = phylogroup, y = count, label = paste0(round(pct_phylogroup, 2), "%")),
    inherit.aes = FALSE,
    vjust = -0.5,
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  geom_text(
    data = ecoli_phylogroups_freq %>%
      filter(pipolin_presence == 0 & phylogroup %in% c("A","E","B2","B1","F","D","C","G")),
    aes(x = phylogroup, y = total, label = paste0(round(pct_phylogroup, 2), "%")),
    inherit.aes = FALSE,
    vjust = -0.8,
    color = "black",
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = seq(0, max(ecoli_phylogroups_freq$count*1.2), 1000)
  ) +
  scale_fill_manual(
    name = "Pipolinas",
    values = c("0" = "#87CEFA", "1" = "#FA8072"), 
    labels = c("Ausencia", "Presencia")) +
  theme_classic() +
  labs(
    title = expression("Distribución de genomas con pipolinas entre filogrupos de " * italic("Escherichia coli")),
    x = "Filogrupo",
    y = "Nº Genomas",
    fill = 'Pipolins'
  ) +
  theme(
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )
plot_phylogroups
#=========================================================================
# Step 2: Plot ST by pipolin presence
#=========================================================================
# Make long table
ecoli_st_freq <- ecoli_dataset %>%
  group_by(ST, pipolin_presence) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ST) %>%
  mutate(
    total = sum(count),
    pct_st = count / total * 100
  ) %>%
  ungroup()

# Plot
plot_ST <- ggplot(ecoli_st_freq %>%
         filter(ST != '-' & total >= 150), 
       aes(x = ST, y = count, fill = factor(pipolin_presence))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    data = ecoli_st_freq %>%
      filter(pipolin_presence == 1 & ST != '-' & total >= 150),
    aes(x = ST, y = count, label = paste0(round(pct_st, 2), "%")),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 2.5,
    color = "black"  
  ) +
  geom_text(
    data = ecoli_st_freq %>%
      filter(pipolin_presence == 0 & ST != '-' & total >= 150),
    aes(x = ST, y = total, label = paste0(round(pct_st, 2), "%")),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 2.5,
    color = "black"
  ) +
  labs(
    x = "ST",
    y = "Num Genomes",
    fill = 'Pipolins'
  ) +
  scale_y_continuous(
    breaks = seq(0, max(ecoli_st_freq$count*1.2), 1000)
  ) +
  scale_fill_manual(
    values = c("0" = "#87CEFA", "1" = "#FA8072"), 
    labels = c("Absence", "Presence")) +
  theme_bw()

#=========================================================================
# Step 3: Plot phylogroups & ST by pipolin presence
#=========================================================================
# Make long table
ecoli_st_phylogroup_freq <- ecoli_dataset %>%
  group_by(ST, phylogroup, pipolin_presence) %>%
  summarise(count = n(), .groups = "drop") %>%
  
  # total genomes per ST
  group_by(ST) %>%
  mutate(
    total_st = sum(count),
    pct_st = count / total_st * 100
  ) %>%
  
  # total genomes per ST in each phylogroup
  group_by(ST, phylogroup) %>%
  mutate(
    total_st_phylogroup = sum(count)
  ) %>%
  
  # total genomes per phylogroup
  group_by(phylogroup) %>%
  mutate(
    total_phylogroup = sum(count)
  ) %>%
  
  # percentage of the phylogroup that belongs to that ST
  mutate(
    pct_st_phylogroup = total_st_phylogroup / total_phylogroup * 100
  ) %>%
  
  ungroup()


# Filtered dataset for the plot
plot_data <- ecoli_st_phylogroup_freq %>%
  filter(
    ST != "-",
    total_st >= 150,
    pct_st_phylogroup >= 1,
    phylogroup %in% c("A","E","B2","B1","F","D","C","G")
  )

# Labels for the facet_wrap
phylo_counts <- ecoli_dataset %>%
  filter(phylogroup %in% c("A","E","B2","B1","F","D","C","G")) %>%
  count(phylogroup, name = "n")

facet_labels <- setNames(
  paste0(phylo_counts$phylogroup, " (n = ", phylo_counts$n, ")"),
  phylo_counts$phylogroup
)

# Plot
plot_phylogrups_ST <- ggplot(
  plot_data,
  aes(
    x = reorder_within(ST, total_st_phylogroup, phylogroup),
    y = count,
    fill = factor(pipolin_presence)
  )
) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    data = plot_data %>% filter(pipolin_presence == 1),
    aes(
      x = reorder_within(ST, total_st_phylogroup, phylogroup),
      y = count,
      label = paste0(round(pct_st, 2), "%")
    ),
    inherit.aes = FALSE,
    vjust = 0.5,
    size = 2.5
  ) +
  geom_text(
    data = plot_data %>% filter(pipolin_presence == 0),
    aes(
      x = reorder_within(ST, total_st_phylogroup, phylogroup),
      y = total_st_phylogroup,
      label = paste0(round(pct_st, 2), "%")
    ),
    inherit.aes = FALSE,
    vjust = -0.3,
    size = 2.5
  ) +
  facet_wrap(~phylogroup, scales = "free_x", labeller = labeller(phylogroup = facet_labels)) +
  scale_x_reordered() +
  labs(
    x = "ST",
    y = "Num Genomes",
    fill = "Pipolins"
  ) +
  scale_y_continuous(
    breaks = seq(0, max(plot_data$count) * 1.5, 500),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(
    values = c("0" = "#87CEFA", "1" = "#FA8072"),
    labels = c("Absence", "Presence")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# ALL PLOTS
plot_phylogroups
plot_ST
plot_phylogrups_ST


#=========================================================================
# Step 4: Extract data
#=========================================================================
### A) Phylogroup A, ST48
A_48 <- ecoli_dataset[ecoli_dataset$ST == "48" & ecoli_dataset$phylogroup == "A",]
#writeLines(A_48$assembly, con = "ecoli_A_48.txt")

### B) Total genomes with phylogroup: file with Assembly - Phylogroup
nrow(ecoli_dataset)
unique(ecoli_dataset$phylogroup)
ecoli_dataset_w_phylogroup <- ecoli_dataset[ecoli_dataset$phylogroup %in% c(
  'A', 'B1', 'B2', 'C', 'D', 'E', 'F', 'G'
),]
unique(ecoli_dataset_w_phylogroup$phylogroup)
nrow(ecoli_dataset_w_phylogroup)

write.table(ecoli_dataset, file = 'ecoli_dataset.txt', sep = '\t', quote = FALSE, row.names = FALSE)
#write.table(ecoli_dataset_w_phylogroup[,c('assembly', 'phylogroup')], file = 'ecoli_genomes_phylogroups.txt', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(ecoli_dataset_w_phylogroup, file = 'ecoli_dataset_phylogroups.txt', sep = '\t', quote = FALSE, row.names = FALSE)

### C) Genomes with pipolins and phylogroup: file with Assembly - Phylogroup
#write.table(ecoli_dataset_w_phylogroup[ecoli_dataset_w_phylogroup$pipolin_presence == TRUE ,c('assembly', 'phylogroup')], file = 'ecoli_pipolin_genomes_phylogroups.txt', sep = '\t', quote = FALSE, row.names = FALSE)
