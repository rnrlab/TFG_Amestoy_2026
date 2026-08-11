#### Graficacion de resultados de ExplorePipolin ####

packages <- c("dplyr", "ggplot2", "shadowtext", "patchwork")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

#=====================================================================================================
#===== Step 0: Importing data ========================================================================
#=====================================================================================================

# IMPORTANT!! --> change to pipolin_metadata_filtered.txt
pipolin_metadata <- read.table(
  file = "pipolin_metadata_filtered.txt",
  sep = "\t",
  header = TRUE
)
pipolin_metadata <- pipolin_metadata %>%
  filter(Completeness >= 90 & Contamination <= 5)

bacteria_summary <- read.table(
  file = "bacteria_summary.txt",
  header = TRUE,
  sep = "\t"
)

#=====================================================================================================
#===== Step 1: Pipolins per assembly + piPolBs per pipolin (figure 2A) ===============================
#=====================================================================================================

## A) Prepare common data
order_more_20 <- pipolin_metadata %>%
  group_by(order) %>%
  summarise(
    total_pipolins = sum(num_of_pipolins, na.rm = TRUE)
  ) %>%
  filter(total_pipolins > 20)
pipolin_order <- pipolin_metadata[pipolin_metadata$order %in% order_more_20$order,]

## B) Pipolins per assembly
# Prepare data
pipolin_order_table <- as.data.frame(table(
  pipolin_order$class,
  pipolin_order$order,
  pipolin_order$num_of_pipolins
))
colnames(pipolin_order_table) <- c("class", "order", "num_of_pipolins", "Freq")

length(unique(pipolin_order$order)) # 14 ordenes diferentes
color_pipolin_orders <- c("#8a2b32", "#fa7223", "#25d9c1", "#5bad51", "#c782ff",
                          "#7FE072", "#8e2abf", "#4781e6", "#eb3f4d","#ff94f1",
                          "#e32bac", "#ffa947", "#1f7a8c", "grey40")
length(unique(pipolin_order$order)) == length(color_pipolin_orders)

class(pipolin_order_table$num_of_pipolins)
pipolin_order_table$num_of_pipolins <- as.integer(as.character(pipolin_order_table$num_of_pipolins))

# Plot
n_o_pipolins <- ggplot(pipolin_order_table[pipolin_order_table$class != "Unknown" & pipolin_order_table$order != "Unknown" & pipolin_order_table$num_of_pipolins>0 & pipolin_order_table$Freq>0,], 
       aes(x = num_of_pipolins, y = Freq + 1, fill = order)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_pipolin_orders, name = "Orden") +
  coord_flip() +
  scale_y_log10(labels = c()) +
  geom_shadowtext(
    aes(label = ifelse(Freq>0, Freq, "")),
    size = 7,
    angle = 0,
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Pipolinas por genoma") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"), axis.ticks.x = element_blank()) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title = element_text(size = 23, face = "bold"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 26)
  )

## C) piPolBs per assembly
# Prepare data
pipolb_order_table <- as.data.frame(table(
  pipolin_order$class,
  pipolin_order$order,
  pipolin_order$num_of_piPolB
))
colnames(pipolb_order_table) <- c("class", "order", "num_of_piPolB", "Freq")

length(unique(pipolin_order$order)) # 14 ordenes diferentes
length(unique(pipolin_order$order)) == length(color_pipolin_orders)

class(pipolb_order_table$num_of_piPolB)
pipolb_order_table$num_of_piPolB <- as.integer(as.character(pipolb_order_table$num_of_piPolB))

# Plot
n_o_pipolbs <- ggplot(pipolb_order_table[pipolb_order_table$class != "Unknown" & pipolb_order_table$order != "Unknown" & pipolb_order_table$num_of_piPolB>0 & pipolb_order_table$Freq>0,], 
                       aes(x = num_of_piPolB, y = Freq + 1, fill = order)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_pipolin_orders, name = "Orden") +
  coord_flip() +
  scale_y_log10(labels = c()) +
  geom_shadowtext(
    aes(label = ifelse(Freq>0, Freq, "")),
    size = 7,
    angle = 0,
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("piPolBs por genoma") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"), axis.ticks.x = element_blank()) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title = element_text(size = 23, face = "bold"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 26)
  )

## D) Join plots and save
(pipolins_pipolbs_per_assembly <- n_o_pipolins / n_o_pipolbs + plot_layout(ncol = 1, nrow = 2, guides = "collect"))

#=====================================================================================================
#===== Step 2: Number of direct repeats in pipolins (figure 2B) ======================================
#=====================================================================================================

#=====================================================================================================
#===== Step 3: Ratio of pipolin reconstruction gaps and number of direct repeats (figure 2C) =========
#=====================================================================================================

#=====================================================================================================
#===== Step 4: Length of pipolins by bacterial order (figure 2D) =====================================
#=====================================================================================================
