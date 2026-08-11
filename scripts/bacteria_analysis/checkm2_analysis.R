#### Analysis of CheckM2 results on 16417 genomes with pipolins ####
"
We will analyze results of CheckM2 output over bacterial genomes with pipolins.
This results include completeness and contamination percentages per genome.

We end by exporting only those genomes with a certain standar level of
completeness and contamination:
 - Completeness >= 90%
 - Contamination <= 5%
"

packages <- c("ggplot2", "dplyr")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}


#=========================================================================
#===== Step 0: Importing data from CheckM2 ===============================
#=========================================================================

checkm2 <- read.table(
  file = "checkm2_all_genomes.tsv",
  sep = "\t",
  header = TRUE
)

#=========================================================================
#===== Step 1: Completeness analysis =====================================
#=========================================================================

## A) Histogram
ggplot(data = checkm2, aes(x = Completeness)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "white") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5)
  ) +
  labs(
    title = "Distribution of genome Completeness",
    x = "Completeness (%)",
    y = "Frecuency"
  )

## B) Histrogram (Completeness < 100)
ggplot(data = checkm2[checkm2$Completeness<90,], aes(x = Completeness)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "white") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5)
  ) +
  labs(
    title = "Distribution of genome Completeness",
    x = "Completeness (%)",
    y = "Frecuency"
  )

#=========================================================================
#===== Step 2: Completeness x Contamination analysis =====================
#=========================================================================

## A) Scatter plot Completeness vs Contamination
ggplot(checkm2, aes(x = Completeness, y = Contamination)) +
  geom_point(color = "skyblue", alpha = 0.6, size = 1.8) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 5)
  ) +
  theme_minimal() +
  labs(
    title = "Genome quality by CheckM2",
    x = "Completeness (%)",
    y = "Contamination (%)"
  )

## B) Double histogram
checkm2_long <- stack(checkm2[,2:3])
colnames(checkm2_long) <- c("value", "type")

checkm2_plot <- ggplot(checkm2_long, aes(x = value, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30, color = "black") +
  scale_fill_manual(
    values = c("steelblue", "firebrick"),
                    
    labels = c("% Completitud", "% Contaminación")
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5)
  ) +
  scale_y_continuous(
    breaks = seq(0, max(nrow(checkm2))*1.5, by = 1000)
  ) +
  labs(
    title = "Validación de la calidad de los genomas por CheckM2",
    x = "Distribución (%)", 
    y = "Número de genomas", 
    fill = "Variable") +
  theme_classic(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 26)
  )
checkm2_plot

#=========================================================================
#===== Step 3: Variable analysis =========================================
#=========================================================================

## A) Completeness summary table
table_completeness <- data.frame(
  Completeness = c("Total", ">=90%", "<90%"),
  Num.Assemblies = c(
    nrow(checkm2),
    nrow(checkm2[checkm2$Completeness>=90,]),
    nrow(checkm2[checkm2$Completeness<90,])
  ),
  Prop.Assemblies = c(
    (nrow(checkm2)/nrow(checkm2))*100,
    (nrow(checkm2[checkm2$Completeness>=90,])/nrow(checkm2))*100,
    (nrow(checkm2[checkm2$Completeness<90,])/nrow(checkm2))*100
  ),
  Min = numeric(3),
  Q1 = numeric(3),
  Median = numeric(3),
  Mean = numeric(3),
  Q3 = numeric(3),
  Max = numeric(3)
)
table_completeness[1, 4:9] <- summary(checkm2$Completeness)
table_completeness[2, 4:9] <-summary(checkm2$Completeness[checkm2$Completeness>=90])
table_completeness[3, 4:9] <-summary(checkm2$Completeness[checkm2$Completeness<90])
table_completeness
write.table(table_completeness, "table_completeness.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## B) Contamination summary table
table_contamination <- data.frame(
  Contamination = c("Total", "<=5%", ">5%"),
  Num.Assemblies = c(
    nrow(checkm2),
    nrow(checkm2[checkm2$Contamination<=5,]),
    nrow(checkm2[checkm2$Contamination>5,])
  ),
  Prop.Assemblies = c(
    (nrow(checkm2)/nrow(checkm2))*100,
    (nrow(checkm2[checkm2$Contamination<=5,])/nrow(checkm2))*100,
    (nrow(checkm2[checkm2$Contamination>5,])/nrow(checkm2))*100
  ),
  Min = numeric(3),
  Q1 = numeric(3),
  Median = numeric(3),
  Mean = numeric(3),
  Q3 = numeric(3),
  Max = numeric(3)
)
table_contamination[1, 4:9] <- summary(checkm2$Contamination)
table_contamination[2, 4:9] <-summary(checkm2$Contamination[checkm2$Contamination<=5])
table_contamination[3, 4:9] <-summary(checkm2$Contamination[checkm2$Contamination>5])
table_contamination
write.table(table_contamination, "table_contamination.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#=========================================================================
#===== Step 4: Export filtered table =====================================
#=========================================================================

checkm2_filtered <- checkm2[checkm2$Completeness>=90 & checkm2$Contamination<=5,]

write.table(
  checkm2_filtered,
  file = "checkm2_filtered_all_genomes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

