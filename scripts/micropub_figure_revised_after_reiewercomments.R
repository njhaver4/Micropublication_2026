
##R version: 4.2.1 library: Rlibs
rm(list = ls())
library(ggplot2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
#install.packages("cowplot")
library(cowplot)
library(RColorBrewer)
#install.packages("patchwork")  
library(patchwork)

data_briggsae<-read.csv(file = "/vast/eande106/projects/Nikita/filtration_optimization_briggsae/Micropublication_github/data/raw_data/briggsae_filtration.csv",header = TRUE,sep = ',')
data_tropicalis<-read.csv(file = "/vast/eande106/projects/Nikita/filtration_optimization_briggsae/Micropublication_github/data/raw_data/tropicalis_filtration.csv",header = TRUE,sep = ',')
data_hta<-read.csv(file = "/vast/eande106/projects/Nikita/filtration_optimization_briggsae/Micropublication_github/data/raw_data/briggsae_HTLDA_easyXpress3.0-processed.csv",header = TRUE,sep = ',')

data_briggsae <- data_briggsae %>%
  mutate(
    Strain = recode(
      Strain,
      "PB420" = "CGC2"
    )
  )

average_embryos_briggsae <- data_briggsae %>%
  group_by(Strain, Condition) %>%
  summarise(
    mean_embryos_per_plate = mean(
      Embryos_per_plate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

data_tropicalis <- data_tropicalis %>%
  mutate(
    Strain = recode(
      Strain,
      "NIC58" = "CGC3"
    )
  )

average_embryos_tropicalis <- data_tropicalis %>%
  group_by(Strain, Experiment) %>%
  summarise(
    mean_embryos_per_plate = mean(
      Embryos_per_plate,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

data_hta <- data_hta %>%
  mutate(
    strain = recode(
      strain,
      "PB420" = "CGC2"
    )
  )

################Significant values###############
##                                             ##
##                                             ##
################Significant values###############

sig_label <- function(p) {
  case_when(
    p <= 0.0001 ~ "****",
    p <= 0.001  ~ "***",
    p <= 0.01   ~ "**",
    p <= 0.05   ~ "*",
    TRUE        ~ "" ##Do not want ns in the plot
  )
}


################HTA figure###############
##                                     ##
##                                     ##
################HTA figure###############
filtered_data_1 <- data_hta %>%
  dplyr::select(Metadata_Experiment, drug, concentration_um, strain, median_wormlength_um)

filtered_data_1 <- filtered_data_1 %>%
  mutate(Metadata_Experiment = recode(Metadata_Experiment,
                                      "CbfiltHTA2" = "A2",
                                      "CbfiltHTA3" = "A3",
                                      "CbfiltHT4" = "A4",
                                      "CbfiltHTA5" = "A5"))


filtered_data_1 <- filtered_data_1 %>%
  mutate(
    condition = sub(".*_", "", drug)  # extract text after underscore
  )

filtered_data_1 <- filtered_data_1 %>%
  mutate(
    person = sub("_.*", "", drug) #extract text before underscore
  )

#####Removing and separating by person
##For MM remove ECA2666 from HTA2
##For MM remove JU3237 from HTA3
##For MM remove NIC1667 from HTA4


filtered_data_2 <- filtered_data_1 %>%
  filter(
    !(Metadata_Experiment == "A2" & strain == "ECA2666" & person == "MM"),
    !(Metadata_Experiment == "A3" & strain == "JU3237" & person == "MM"),
    !(Metadata_Experiment == "A4" & strain == "NIC1667" & person == "MM")
  )


####Wilcoxon tests per strain (ignoring person)
wilcox_results_1 <- filtered_data_2 %>%
  group_by(strain) %>%
  summarise(
    p_value = wilcox.test(
      median_wormlength_um ~ condition
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sig_label(p_value)
  )

####Label positions: one per strain, dynamic x-position
label_positions_1 <- filtered_data_2 %>%
  group_by(strain) %>%
  summarise(
    x_pos = 1.5,
    y_pos = 830,
    .groups = "drop"
  ) %>%
  left_join(
    wilcox_results_1,
    by = "strain"
  )

#plot
p_hta <- ggplot(
  filtered_data_2,
  aes(
    x = condition,
    y = median_wormlength_um
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(
      group = interaction(
        strain, condition
      )
    ),
    fill = "grey90",
    color = "black",
    width = 0.55,
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  
  # Individual observations
  geom_jitter(
    aes(color = person),
    width = 0.12,
    alpha = 0.6,
    size = 0.6
  ) +
  
  scale_color_manual(
    values = c(
      "#0072B2",  # blue
      "#E69F00"   # orange
    )
  ) +
  
  # Significance labels
  geom_text(
    data = label_positions_1,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    family = "Helvetica",
    fontface = "italic"
  ) +
  
  # Facet by strain
  facet_grid(
    . ~ strain,
    scales = "fixed"
  ) +
  
  # X-axis labels
  scale_x_discrete(
    labels = c(
      "3% PEG" = "PEG",
      "3% PEG with 0.05% T" = "PEG + Tween 20"
    )
  ) +
  
  # Give annotations room
  scale_y_continuous(
    breaks = seq(400, 800, by = 100),
    expand = expansion(
      mult = c(0.02, 0.05)
    ),
    sec.axis = sec_axis(
      ~ .,
      name = expression(italic("C. briggsae")),
      breaks = NULL
    )
  ) +
  
  coord_cartesian(
    ylim = c(400, 850)
  ) +
  
  labs(
    x = "Condition",
    y = "Median animal length (µm)",
    color = "Experimentalist"
  ) +
  
  theme_classic(
    base_family = "Helvetica"
  ) +
  
  theme(
    
    axis.title.x = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    # UNBOLD y-axis title
    axis.title.y = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    axis.text.x = element_text(
      family = "Helvetica",
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      family = "Helvetica",
      size = 10,
      color = "black"
    ),
    
    # UNBOLD strain labels
    strip.text = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.title = element_text(
      family = "Helvetica",
      size = 10,
      face = "plain"
    ),
    
    legend.text = element_text(
      family = "Helvetica",
      size = 10
    ),
    
    legend.position = "none",
    
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      0.5,
      0.5,
      0.5,
      0.5
    )
  )

p_hta
################C. briggsae###############
##                                      ##
##                                      ##
################C. briggsae###############

wilcox_results_2 <- data_briggsae %>%
  group_by(Strain) %>%
  summarise(
    p_value = wilcox.test(
      Embryos_per_plate ~ Condition
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sig_label(p_value)
  )

label_positions_2 <- data_briggsae %>%
  group_by(Strain) %>%
  summarise(
    y_pos = 10000,
    x_pos = mean(as.numeric(factor(Condition))),
    .groups = "drop"
  ) %>%
  left_join(wilcox_results_2, by = "Strain")


p_briggsae <- ggplot(
  data_briggsae,
  aes(
    x = Condition,
    y = Embryos_per_plate
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(
      group = interaction(
        Strain,
        Condition
      )
    ),
    fill = "grey90",
    color = "black",
    width = 0.55,
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  
  # Individual observations
  geom_jitter(
    aes(color = Person),
    width = 0.12,
    alpha = 0.6,
    size = 1.0
  ) +
  
  scale_color_manual(
    values = c(
      "#0072B2",  # blue
      "#E69F00"   # orange
    )
  ) +
  # Significance labels
  geom_text(
    data = label_positions_2,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    family = "Helvetica",
    fontface = "italic"
  ) +
  
  # Facet by strain
  facet_grid(
    . ~ Strain,
    scales = "fixed"
  ) +
  
  # X-axis labels
  scale_x_discrete(
    labels = c(
      "3% PEG" = "PEG",
      "3% PEG with 0.05% Tween" = "PEG + Tween 20"
    )
  ) +
  
  # Give annotations room
  scale_y_continuous(
    breaks = seq(0, 10000, by = 2000),
    expand = expansion(
      mult = c(0.02, 0.20)
    ),
    sec.axis = sec_axis(
      ~ .,
      name = expression(italic("C. briggsae")),
      breaks = NULL
    )
  ) +
  
  #scale_y_continuous(
   # breaks = seq(0, 10000, by = 2000),
    #expand = expansion(
    #  mult = c(0.02, 0.20)
   # )
 # ) +
  
  labs(
    x = "Condition",
    y = "Embryos per 6 cm plate",
    color = "Experimentalist"
  ) +
  
  theme_classic(
    base_family = "Helvetica"
  ) +
  
  theme(
    
    axis.title.x = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    # UNBOLD y-axis title
    axis.title.y = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    axis.text.x = element_text(
      family = "Helvetica",
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      family = "Helvetica",
      size = 10,
      color = "black"
    ),
    
    # UNBOLD strain labels
    strip.text = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.title = element_text(
      family = "Helvetica",
      size = 10,
      face = "plain"
    ),
    
    legend.text = element_text(
      family = "Helvetica",
      size = 10
    ),
    
    legend.position = "none",
    
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      0.5,
      0.5,
      0.5,
      0.5
    )
  )

p_briggsae

###############C. tropicalis###############
##                                        ##
##                                        ##
################C. tropicalis###############

wilcox_results_3 <- data_tropicalis %>%
  group_by(Strain) %>%
  summarise(
    p_value = wilcox.test(
      Embryos_per_plate ~ Experiment
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sig_label(p_value)
  )


label_positions_3 <- data_tropicalis %>%
  group_by(Strain) %>%
  summarise(
    y_pos = 9000,
    x_pos = mean(as.numeric(factor(Experiment))),
    .groups = "drop"
  ) %>%
  left_join(wilcox_results_3, by = "Strain")   



p_tropicalis <- ggplot(
  data_tropicalis,
  aes(
    x = Experiment,
    y = Embryos_per_plate
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(
      group = interaction(
        Strain,
        Experiment
      )
    ),
    fill = "grey90",
    color = "black",
    width = 0.55,
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  
  # Individual observations
  geom_jitter(
    aes(color = Person),
    width = 0.12,
    alpha = 0.6,
    size = 1.0
  ) +
  
  scale_color_manual(
    values = c(
      "#0072B2",  # blue
      "#E69F00"   # orange
    )
  ) +
  # Significance labels
  geom_text(
    data = label_positions_3,
    aes(
      x = x_pos,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    family = "Helvetica",
    fontface = "italic"
  ) +
  
  # Facet by strain
  facet_grid(
    . ~ Strain,
    scales = "free"
  ) +
  
  # X-axis labels
  scale_x_discrete(
    labels = c(
      "PEG" = "PEG",
      "PEG_Twen20" = "PEG + Tween 20"
    )
  ) +
  # Give annotations room
  scale_y_continuous(
    breaks = seq(0, 10000, by = 2000),
    expand = expansion(
      mult = c(0.02, 0.20)
    ),
    sec.axis = sec_axis(
      ~ .,
      name = expression(italic("C. tropicalis")),
      breaks = NULL
    )
  ) + 
  
  labs(
    x = "Condition",
    y = "Embryos per 6 cm plate",
    color = "Experimentalist"
  ) +
  
  theme_classic(
    base_family = "Helvetica"
  ) +
  
  theme(
    
    axis.title.x = element_text(
      family = "Helvetica",
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    # UNBOLD y-axis title
    axis.title.y = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    axis.text.x = element_text(
      family = "Helvetica",
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      family = "Helvetica",
      size = 10,
      color = "black"
    ),
    
    # UNBOLD strain labels
    strip.text = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.title = element_text(
      family = "Helvetica",
      size = 10,
      face = "plain"
    ),
    
    legend.text = element_text(
      family = "Helvetica",
      size = 10
    ),
    
    legend.position = "none",
    
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      0.5,
      0.5,
      0.5,
      0.5
    )
  )

p_tropicalis

################final_figure###############
##                                       ##
##                                       ##
################final_figure###############

########## Panel A ##########
p_briggsae_final <- p_briggsae +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "none",
    
    plot.margin = margin(
      0.5, 0.5, 0.5, 0.5
    )
  )

p_briggsae_final

########## Panel B ##########
p_hta_final <- p_hta +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "none",
    
    plot.margin = margin(
      0.5, 0.5, 0.5, 0.5
    )
  )

p_hta_final

########## Panel C: C. tropicalis ##########
p_tropicalis_final <- p_tropicalis +
  theme(
    axis.title.x = element_blank(), #Remove the "Condition" label
    
    # KEEP x-axis labels
    axis.text.x = element_text(
      family = "Helvetica",
      size = 10,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    # KEEP x-axis ticks
    axis.ticks.x = element_line(
      color = "black"
    ),
    
    axis.title.y = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    strip.text = element_text(
      family = "Helvetica",
      size = 12,
      face = "plain",
      color = "black"
    ),
    
    # REMOVE legend
    legend.position = "none",
    
    plot.margin = margin(
      0.5, 0.5, 0.5, 0.5
    )
  )


p_tropicalis_final

########## Combined figure ##########
library(patchwork)
combined_figure <-
  (p_briggsae_final /
     p_hta_final /
     p_tropicalis_final) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      face = "bold",
      size = 14,
      color = "black"
    )
  )

combined_figure


ggsave("/vast/eande106/projects/Nikita/filtration_optimization_briggsae/Micropublication_github/plots/20261103.png", 
plot = combined_figure, width = 8, height = 11, unit = "in", dpi = 600)

