#### Analysis of the Prevalence of Bacterial Genomes with Pipolins ####

"
We will use the file bacteria_summary.txt as input. 
This file contains the bacterial genera with genomes containing pipolins.
It includes their Family, Order, and Class.
It also includes the total number of genomes with pipolins and their prevalence.
"

packages <- c("dplyr", "ggplot2", "forcats")
for(package in packages){
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

c <- read.table("datasets_for_benchling/pipolin_metadata_filtered6.txt", sep = "\t", header = T)
c2 <- c[c$Completeness >= 90 & c$Contamination <= 5,]
summary(c2$Completeness[c2$species == "Escherichia_coli"])

#=========================================================================
#===== Step 0: Importing and Processing data =============================
#=========================================================================

bacteria_summary <- read.table(
  file = "bacteria_summary.txt",
  header = TRUE,
  sep = "\t"
)

bacteria_summary <- bacteria_summary %>%
  arrange(Class,Order,Family,Genus)

bacteria_summary$Pipol.Perc <- bacteria_summary$Pipol.Prop * 100

#=========================================================================
#===== Step 1: Prevalence Plot ===========================================
#=========================================================================
nrow(bacteria_summary[bacteria_summary$Pipol.Genomes>3,])
nrow(bacteria_summary[bacteria_summary$Pipol.Genomes>4,])

install.packages("ggtext")
library(ggtext)

ggplot(bacteria_summary %>%
         filter(Pipol.Genomes>4 & Pipol.Prevalence > 1), 
       aes(x = Genus, y = Pipol.Prevalence, group = Order)) +
  geom_bar(aes(fill = Order), stat = "identity", col = "black", width = 1) +
  scale_fill_discrete(name = "Orden") +
  scale_y_continuous(expand = c(0,0), limits=c(0,100), breaks=c(20,40,60,80,100)) +
  geom_richtext(
    aes(
      label = ifelse(
        Pipol.Prevalence > 50,
        # Etiquetas DENTRO de la barra
        paste0(
          "<span style='color:white;'>", round(Pipol.Prevalence,2), "%</span>",
          "<span style='color:black;'>  (", Total.Genomes, ")</span>"
        ),
        # Etiquetas FUERA de la barra
        paste0(
          "<span style='color:black;'>", round(Pipol.Prevalence,2), "%</span>",
          "<span style='color:#777777;'>  (", Total.Genomes, ")</span>"
        )
      ),
      hjust = ifelse(Pipol.Prevalence > 50, 1.05, -0.04)  # dentro centrado, fuera a la derecha
    ),
    angle = 90,
    size = 4.5,
    fontface = "bold",
    fill = NA,          # sin fondo
    label.color = NA,   # sin borde
    vjust = 0.5
  ) +
  facet_grid(.~fct_reorder(Order,Class), scale = "free_x", space = "free") +
  theme_linedraw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1.03, vjust = 1.05, face = "italic", size = 14)) +
  theme(
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.y = element_text(size = 15),
    title = element_text(size = 20),
    legend.text = element_text(size = 15)
  ) +
  theme(plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0.4, unit = "cm"), legend.text = element_text(face = "italic"), legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.text.x.top = element_blank(), panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white")) +
  labs(title = "Distribución de pipolinas por género bacteriano", x = "Género", y = "Prevalencia %") +
  coord_cartesian(clip = "off")
