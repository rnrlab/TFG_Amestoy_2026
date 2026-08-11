#### Analysis of CONJScan results over 16417 genomes with pipolins ####

"
We analyzed the CONJScan results (MacSyFinder). 
This includes both conjugative element profiles (HMM) 
and plasmid profiles (MOB relaxases).
"

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "purrr", "colorspace", "patchwork", "tidytext")
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
# Assemblys to CONJScan systems detected:
systems_per_genome <- read.table(
  file = 'systems_per_genome.txt',
  sep = '\t',
  header = TRUE
)

# Assemblies with pipolins + info
pipolin_metadata_filtered <- read.table(
  file = "pipolin_metadata_filtered.txt",
  sep = '\t',
  header = TRUE
)

# MOB types presence-absence matrix
mobs_matrix <- read.table(
  file = 'mobs_matrix.txt',
  sep = '\t',
  header = TRUE
)

# MOBs per genome
mobs_per_genome <- read.table(
  file = 'mobs_per_genome.txt',
  sep = '\t',
  header = TRUE
)

## B) Process
# Remove extra or redundant columns [replicon]
systems_per_genome$replicon <- NULL

# Merge tables by assembly
pipolin_metadata_filtered <- merge(pipolin_metadata_filtered, systems_per_genome, by = 'assembly')
pipolin_metadata_filtered <- merge(pipolin_metadata_filtered, mobs_matrix, by = 'assembly')
pipolin_metadata_filtered <- merge(pipolin_metadata_filtered, mobs_per_genome, by = 'assembly')

# Filter by CheckM2 standars
sum(pipolin_metadata_filtered$Completeness >= 90 & pipolin_metadata_filtered$Contamination <= 5)
pipolin_metadata_filtered_original <- pipolin_metadata_filtered
pipolin_metadata_filtered <- pipolin_metadata_filtered[pipolin_metadata_filtered$Completeness >= 90 & pipolin_metadata_filtered$Contamination <= 5, ]

#=========================================================================
#===== Step 1: Distribution of all systems ===============================
#=========================================================================

### A) Frequency of all systems in all genomes
all_systems <- data.frame(
  system = names(colSums(pipolin_metadata_filtered[34:67])),
  Freq = colSums(pipolin_metadata_filtered[34:67])
)
rownames(all_systems) <- NULL

# Plot only Freqs
ggplot(data = all_systems[order(all_systems$Freq, decreasing = TRUE),], aes(x = reorder(system, Freq), y = Freq)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme_minimal()

### B) Frequency of all systems in all genomes by genus 
# Combine Chromosome and Plasmid
all_systems_genus <- pipolin_metadata_filtered %>%
  pivot_longer(
    cols = 34:67,
    names_to = 'system',
    values_to = 'value'
  ) %>%
  mutate(system = gsub("CONJScan\\.(Plasmids|Chromosome)\\.", "CONJScan.", system)) %>%
  group_by(class, order, genus, system) %>%
  summarise(total_systems = sum(value, na.rm = TRUE), .groups = 'drop')

# Filter for plot
all_systems_genus_plot <- all_systems_genus[all_systems_genus$total_systems > 10,]

# Organize
all_systems_genus_plot <- all_systems_genus_plot %>%
  mutate(system = sub(".*\\.", "", system))

all_systems_genus_plot <- all_systems_genus_plot %>%
  arrange(order, genus) %>%
  mutate(genus = factor(genus, levels = unique(genus)))

# Plot by genus
ggplot(
  data = all_systems_genus_plot %>% 
    filter(total_systems > 100),
  aes(x = reorder(system, total_systems), 
      y = total_systems, 
      fill = genus)
) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  facet_wrap(. ~ genus, scales = "free") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Aeromonas" = "#FA8072",
    "Staphylococcus" = "#ffa947",
    "Citrobacter" = "#5bad51",
    "Enterobacter" = "#7FE072",
    "Escherichia" = "#25d9c1",
    "Limosilactobacillus" = "#4781e6",
    "Pseudosulfitobacter" = "#c782ff",
    "Vibrio" = "#e32bac"
  )) +
  labs(
    title = "Distribución de elementos conjugativos por género",
    x = "Sistema conjugativo",
    y = "Número de sistemas",
    fill = "Género"
  ) +
  theme_classic() +
  theme(
    legend.position = 'none',
    legend.text = element_text(face = "italic", size = 18),
    strip.text = element_text(size = 22, face = "italic"),
    title = element_text(size = 25),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
  )

### C) Percentage of all systems in all genomes by genus 
## Número total de genomas por género
genus_totals <- pipolin_metadata_filtered %>%
  group_by(genus) %>%
  summarise(
    total_genomes = n(),
    .groups = "drop"
  )

## Presencia de sistemas por género
all_systems_genus_perc <- pipolin_metadata_filtered %>%
  
  # Pasar de formato ancho a largo
  pivot_longer(
    cols = 34:67,
    names_to = "system",
    values_to = "value"
  ) %>%
  
  # Unificar nombres de sistemas
  mutate(
    system = gsub(
      "CONJScan\\.(Plasmids|Chromosome)\\.",
      "CONJScan.",
      system
    )
  ) %>%
  
  # Convertir a presencia/ausencia
  mutate(
    value = ifelse(value > 0, 1, 0)
  ) %>%
  
  # EVITAR DUPLICADOS:
  # un mismo genoma solo puede contar una vez
  # para cada sistema
  group_by(assembly, class, order, genus, system) %>%
  summarise(
    value = max(value),
    .groups = "drop"
  ) %>%
  
  # Contar genomas con cada sistema
  group_by(class, order, genus, system) %>%
  summarise(
    genomes_with_system = sum(value),
    .groups = "drop"
  ) %>%
  
  # Añadir total de genomas por género
  left_join(
    genus_totals,
    by = "genus"
  ) %>%
  
  # Calcular porcentaje
  mutate(
    percentage = (genomes_with_system / total_genomes) * 100
  ) %>%
  
  # Limpiar nombres de sistemas
  mutate(
    system = sub(".*\\.", "", system)
  )

all_systems_genus_perc %>%
  filter(percentage > 100)

## Plot by genus
ggplot(
  data = all_systems_genus_perc %>%
    filter(genomes_with_system > 100) %>%
    filter(percentage > 10),
  aes(
    x = reorder_within(system, percentage, genus),
    y = percentage,
    fill = genus
  )
) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(. ~ genus, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  
  scale_fill_manual(values = c(
    "Aeromonas" = "#FA8072",
    "Staphylococcus" = "#ffa947",
    "Citrobacter" = "#5bad51",
    "Enterobacter" = "#7FE072",
    "Escherichia" = "#25d9c1",
    "Limosilactobacillus" = "#4781e6",
    "Pseudosulfitobacter" = "#c782ff",
    "Vibrio" = "#e32bac"
  )) +
  
  labs(
    title = "Porcentaje de genomas con sistemas conjugativos por género",
    x = "Sistema conjugativo",
    y = "Genomas con el sistema (%)",
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
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18)
  )

#=========================================================================
#===== Step 2: Distribution of MOB systems ===============================
#=========================================================================
### A) Frequency of each MOB system
# Define all possible MOB profiles
mobs_profiles <- c('T4SS_MOBB', 'T4SS_MOBC', 'T4SS_MOBF', 'T4SS_MOBH', 'T4SS_MOBM', 'T4SS_MOBP1', 'T4SS_MOBP2', 'T4SS_MOBP3', 'T4SS_MOBQ', 'T4SS_MOBT', 'T4SS_MOBV')

# Add profiles w/o hits
missing_mobs <- setdiff(mobs_profiles, colnames(pipolin_metadata_filtered))
pipolin_metadata_filtered[missing_mobs] <- 0

# Extract MOB presence-absence info from principal dataset
mobs_dataset <- pipolin_metadata_filtered[,c(
  'assembly',
  'organism_name',
  'species',
  'genus',
  'family',
  'order',
  'class',
  mobs_profiles
)]

mobs_dataset_freq <- data.frame(
  mob_system = names(colSums(mobs_dataset[mobs_profiles])),
  Freq = colSums(mobs_dataset[mobs_profiles])
)
rownames(mobs_dataset_freq) <- NULL

# Plot only Freqs
ggplot(data = mobs_dataset_freq[order(mobs_dataset_freq$Freq, decreasing = TRUE),], aes(x = reorder(mob_system, Freq), y = Freq)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme_minimal()

### B) Frequency of each MOB system by genus and order
mobs_dataset_genus <- mobs_dataset %>%
  pivot_longer(
    cols = all_of(mobs_profiles),
    names_to = 'mob_system',
    values_to = 'value'
  ) %>%
  group_by(class, order, genus, mob_system) %>%
  summarise(total_systems = sum(value), .groups = 'drop')

mobs_dataset_genus[c('class', 'order', 'genus', 'mob_system')] <- lapply(mobs_dataset_genus[c('class', 'order', 'genus', 'mob_system')], as.factor)

# Calculate total per MOB system (just to order from more to less frequent)
mob_totals <- mobs_dataset_genus %>%
  group_by(mob_system) %>%
  summarise(total = sum(total_systems), .groups = 'drop') %>%
  filter(total > 5) %>%
  arrange(desc(-total))

# Plot by genus
ggplot(data = mobs_dataset_genus %>%
         filter(total_systems > 5) %>% # Only MOBs with hits in more than 5 genomes (for better visualization)
         mutate(mob_system = factor(mob_system, levels = mob_totals$mob_system)), 
       aes(x = mob_system, y = total_systems, fill = genus)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'right', axis.text.y = element_text(size = 7))

# Plot by genus and order
ggplot(data = mobs_dataset_genus %>%
         filter(total_systems > 30) %>%
         mutate(mob_system = factor(mob_system, levels = mob_totals$mob_system)), 
       aes(x = mob_system, y = total_systems, fill = genus)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  facet_wrap(.~order, scales = "free") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = 'right', axis.text.y = element_text(size = 7)) +
  ggtitle("Distribución de relaxasas MOB por género y orden bacterianos") +
  xlab("Relaxasa MOB") +
  ylab("Número de sistemas") +
  labs(fill = "Género") +
  theme(
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 20),
    legend.text = element_text(size = 15, face = "italic")
  )

### C) Percentage of each MOB system by genus and order
## Total de genomas
genus_totals_mob <- mobs_dataset %>%
  group_by(genus) %>%
  summarise(total_genomes = n(), .groups = "drop")

## Presencia de MOBs por género
mobs_dataset_perc <- mobs_dataset %>%
  pivot_longer(
    cols = all_of(mobs_profiles),
    names_to = "mob_system",
    values_to = "value"
  ) %>%
  
  # Presencia ausencia
  mutate(value = ifelse(value > 0, 1, 0)) %>%
  
  group_by(class, order, genus, mob_system) %>%
  summarise(
    genomes_with_mob = sum(value),
    .groups = "drop"
  ) %>%
  
  # Añadir genomas totales por género y calcular %
  left_join(genus_totals, by = "genus") %>%
  
  mutate(
    percentage = (genomes_with_mob / total_genomes) * 100
  )

## Recalcular orden
mob_totals <- mobs_dataset_perc %>%
  group_by(mob_system) %>%
  summarise(total = mean(percentage), .groups = "drop") %>%
  filter(total > 1) %>%
  arrange(desc(total))

## Plot
ggplot(
  data = mobs_dataset_perc %>%
    filter(genomes_with_mob > 50) %>%
    filter(percentage > 10) %>%
    mutate(mob_system = factor(mob_system, levels = mob_totals$mob_system)),
  aes(
    x = reorder_within(mob_system, percentage, genus),
    y = percentage, 
    fill = genus)
  ) +
  
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black"
  ) +
  
  facet_wrap(. ~ genus, scales = "free") +
  
  coord_flip() +
  
  scale_x_reordered() +
  
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
    "Salmonella" = "#ff94f1"
  )) +
  
  
  labs(
    title = "Porcentaje de genomas con relaxasas MOB por género",
    x = "Relaxasa MOB",
    y = "Genomas con el sistema (%)",
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
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 18)
  )


#=========================================================================
#===== Step 4: Plots only for Escherichia coli ===========================
#=========================================================================
### A) Create dataset for E. coli
pipolin_ecoli <- pipolin_metadata_filtered %>%
  filter(species == "Escherichia_coli")


### B) CONJ systems plot
# Frequencies
all_systems_ecoli <- pipolin_ecoli %>%
  pivot_longer(
    cols = 34:67,
    names_to = 'system',
    values_to = 'value'
  ) %>%
  
  # 1. Unificar Chromosome/Plasmid
  mutate(system = gsub("CONJScan\\.(Plasmids|Chromosome)\\.", "CONJScan.", system)) %>%
  
  # 2. Limpiar nombre
  mutate(system = sub(".*\\.", "", system)) %>%
  
  # 3. Sumar ocurrencias (lo que tú quieres)
  group_by(system) %>%
  summarise(total_systems = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  
  arrange(desc(total_systems))

# Plot
n_genomes <- nrow(pipolin_ecoli)
all_systems_ecoli$Freq_rel <- all_systems_ecoli$total_systems / n_genomes

plot_conj <- ggplot(all_systems_ecoli %>% filter(Freq_rel>0), 
                    aes(x = reorder(system, Freq_rel), 
                        y = Freq_rel)) +
  geom_bar(stat = 'identity', fill = "#4C72B0", color = "black") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Sistemas conjugativos en E. coli",
    x = "Sistema (T4SS / dCONJ)",
    y = "Frecuencia relativa"
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 16),
  )

### C) MOB systems
mobs_profiles <- c(
  'T4SS_MOBB','T4SS_MOBC','T4SS_MOBF','T4SS_MOBH',
  'T4SS_MOBM','T4SS_MOBP1','T4SS_MOBP2','T4SS_MOBP3',
  'T4SS_MOBQ','T4SS_MOBT','T4SS_MOBV'
)

# Añadir columnas faltantes
missing_mobs <- setdiff(mobs_profiles, colnames(pipolin_ecoli))
pipolin_ecoli[missing_mobs] <- 0

mobs_ecoli <- pipolin_ecoli %>%
  pivot_longer(
    cols = all_of(mobs_profiles),
    names_to = 'mob',
    values_to = 'value'
  ) %>%
  
  # limpiar nombre
  mutate(mob = sub("T4SS_", "", mob)) %>%
  
  group_by(mob) %>%
  summarise(total_systems = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  
  arrange(desc(total_systems))

# Plot
mobs_ecoli$Freq_rel <- mobs_ecoli$total_systems / n_genomes

plot_mob <- ggplot(mobs_ecoli %>% filter(Freq_rel>0), 
                   aes(x = reorder(mob, Freq_rel), 
                       y = Freq_rel)) +
  geom_bar(stat = 'identity', fill = "#DD8452", color = "black") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Relaxasas MOB en E. coli",
    x = "Familia MOB",
    y = "Frecuencia relativa"
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    title = element_text(size = 16),
  )

### D) Combine plots
combined_plot <- plot_conj / plot_mob + 
  plot_layout(heights = c(1, 1))

combined_plot
#=========================================================================
#===== Step 3: Exporting data ============================================
#=========================================================================

# Clean dataset
pipolin_metadata_filtered <- pipolin_metadata_filtered %>%
  relocate(T4SS_MOBP3, .after = T4SS_MOBP2)

# Save big dataset to file
write.table(
  pipolin_metadata_filtered,
  file = 'pipolin_metadata_filtered2.txt',
  sep = '\t',
  row.names = FALSE
)
