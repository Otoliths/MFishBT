################################################################################################################################
#                                                                                                                              #
# Purpose:       Fig.2 Plotting                                                                                                #
#                                                                                                                              #
# Author:        Liuyong Ding et al.                                                                                           #
# Contact:       ly_ding@126.com                                                                                               #
# Client:        Liuyong Ding et al.                                                                                           #
#                                                                                                                              #
# Code created:  2022-07-06                                                                                                    #
# Last updated:  2022-07-06                                                                                                    #
# Source:        FishBT_2022                                                                                                   #
#                                                                                                                              #
# Comment:       Sankey diagram representing the proportion of migratory fish studies using biogeochemical tags                #
#                                                                                                                              #
################################################################################################################################

############################### Packages loading ###############################
rm(list = ls())
options(warn = -1)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(ggpubr)

################################# Data loading #################################
dt <- readRDS("input/sankey.rds")

################################ Data cleaning #################################
dt$IUCN <- ifelse(dt$IUCN == "NA", "Not available", dt$IUCN)
names(dt) <- c(
  "Pdf_id",
  "Binomial_fishbase",
  "Biogeochemical tag",
  "Area",
  "GROMS category",
  "IUCN threatened category",
  "Life history stage",
  "Tissues",
  "Scan mode",
  "Abiotic reference",
  "Commercial importance"
)

dt$Tissues <- ifelse(dt$Tissues == "Vertebrae", "Vertebra", dt$Tissues)
dt$Tissues <- ifelse(dt$Tissues == "Teeth", "Tooth", dt$Tissues)
dt$Tissues <- ifelse(dt$Tissues == "Eye lense", "Eye lens", dt$Tissues)
dt$`Life history stage` <- ifelse(dt$`Life history stage` == "Larvae", "Larva", dt$`Life history stage`)
dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "LR", "NT", dt$`IUCN threatened category`)

table(dt$`Biogeochemical tag`) / nrow(dt)
table(dt$Area) / nrow(dt)
table(dt$`GROMS category`) / nrow(dt)
table(dt$`IUCN threatened category`) / nrow(dt)
table(dt$`Commercial importance`) / nrow(dt)
table(dt$`Life history stage`) / nrow(dt)
table(dt$Tissues) / nrow(dt)
table(dt$`Scan mode`) / nrow(dt)
table(dt$`Abiotic reference`) / nrow(dt)

dt$`Biogeochemical tag` <- ifelse(dt$`Biogeochemical tag` == "Microchemistry", "Microchemistry\n64%", dt$`Biogeochemical tag`)
dt$`Biogeochemical tag` <- ifelse(dt$`Biogeochemical tag` == "Isotope", "Isotope\n36%", dt$`Biogeochemical tag`)
dt$Area <- ifelse(dt$Area == "Marine", "Marine\n44%", dt$Area)
dt$Area <- ifelse(dt$Area == "Freshwater", "Inland\n56%", dt$Area)
dt$`GROMS category` <- ifelse(dt$`GROMS category` == "Anadromous", "Anadromous\n33%", dt$`GROMS category`)
dt$`GROMS category` <- ifelse(dt$`GROMS category` == "Oceanodromous", "Oceanodromous\n32%", dt$`GROMS category`)
dt$`GROMS category` <- ifelse(dt$`GROMS category` == "Catadromous", "Catadromous\n17%", dt$`GROMS category`)

dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "LC", "LC\n48%", dt$`IUCN threatened category`)
dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "Not available", "Not available\n21%", dt$`IUCN threatened category`)
dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "EN", "EN\n8%", dt$`IUCN threatened category`)
dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "VU", "VU\n7%", dt$`IUCN threatened category`)
dt$`IUCN threatened category` <- ifelse(dt$`IUCN threatened category` == "CR", "CR\n5%", dt$`IUCN threatened category`)

dt$`Commercial importance` <- ifelse(dt$`Commercial importance` == "Highly commercial", "Highly commercial\n38%", dt$`Commercial importance`)
dt$`Commercial importance` <- ifelse(dt$`Commercial importance` == "Commercial", "Commercial\n43%", dt$`Commercial importance`)
dt$`Commercial importance` <- ifelse(dt$`Commercial importance` == "Minor commercial", "Minor commercial\n11%", dt$`Commercial importance`)


dt$`Life history stage` <- ifelse(dt$`Life history stage` == "Mixed", "Mixed\n48%", dt$`Life history stage`)
dt$`Life history stage` <- ifelse(dt$`Life history stage` == "Adult", "Adult\n24%", dt$`Life history stage`)
dt$`Life history stage` <- ifelse(dt$`Life history stage` == "Juvenile", "Juvenile\n24%", dt$`Life history stage`)

dt$Tissues <- ifelse(dt$Tissues == "Otolith", "Otolith\n77%", dt$Tissues)
dt$Tissues <- ifelse(dt$Tissues == "Muscle", "Muscle\n12%", dt$Tissues)
dt$Tissues <- ifelse(dt$Tissues == "Scale", "Scale\n4%", dt$Tissues)

dt$`Scan mode` <- ifelse(dt$`Scan mode` == "Point", "Point\n62%", dt$`Scan mode`)
dt$`Scan mode` <- ifelse(dt$`Scan mode` == "Line", "Line\n33%", dt$`Scan mode`)
dt$`Scan mode` <- ifelse(dt$`Scan mode` == "Map", "Map\n5%", dt$`Scan mode`)

dt$`Abiotic reference` <- ifelse(dt$`Abiotic reference` == "YES", "YES\n29%", dt$`Abiotic reference`)
dt$`Abiotic reference` <- ifelse(dt$`Abiotic reference` == "NO", "NO\n71%", dt$`Abiotic reference`)


################################### Ploting ####################################
dt %>%
  distinct() %>%
  make_long(`Biogeochemical tag`, Area, `GROMS category`, `IUCN threatened category`, `Commercial importance`, `Life history stage`, Tissues, `Scan mode`, `Abiotic reference`) %>%
  mutate(node = factor(node, levels = c(
    "Isotope\n36%", "Microchemistry\n64%",
    "Marine\n44%", "Inland\n56%",
    "Amphidromous", "Anadromous\n33%", "Catadromous\n17%", "Oceanodromous\n32%", "Potamodromous",
    "Not available\n21%", "DD", "LC\n48%", "NT", "VU\n7%", "EN\n8%", "CR\n5%", "EW", "EX",
    "Unknown", "Of no interest", "Subsistence fisheries", "Minor commercial\n11%", "Commercial\n43%", "Highly commercial\n38%",
    "Mixed\n48%", "Adult\n24%", "Juvenile\n24%", "Larva", "Egg",
    "Others", "Eye lens", "Muscle\n12%", "Tooth", "Vertebra", "Otolith\n77%", "Scale\n4%", "Fin ray",
    "Map\n5%", "Line\n33%", "Point\n62%",
    "NO\n71%", "YES\n29%"
  ))) %>%
  ggplot(aes(
    x = x, next_x = next_x,
    node = node, next_node = next_node,
    label = node
  )) + # label = paste0(node,"\n",percent)
  geom_alluvial(aes(fill = x),
    flow.alpha = .5,
    width = 0.7,
    smooth = 3,
    space = 30,
    na.rm = T
  ) +
  geom_alluvial_text(size = 1.6, color = "black", space = 30) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#B15928", "#FDBF6F", "#FF7F00", "#CAB2D6")) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 2150)) + # label = scales::percent_format(scale = 100)
  scale_x_discrete(expand = c(0.06, 0.01)) +
  labs(x = NULL, y = "Number of studies") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks = element_line(linetype = "blank"),
    axis.text.x = element_text(colour = "black", face = "bold", size = 5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(b = 0.2, unit = "cm")
  )

################################# Fig.2 saving #################################
ggsave("output/figure2.pdf", dpi = 600, width = 20, height = 12, units = "cm")
