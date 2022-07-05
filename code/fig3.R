################################################################################################################################
#                                                                                                                              #
# Purpose:       Fig.3 Plotting                                                                                                #
#                                                                                                                              #
# Author:        Liuyong Ding et al.                                                                                           #
# Contact:       ly_ding@126.com                                                                                               #
# Client:        Liuyong Ding et al.                                                                                           #
#                                                                                                                              #
# Code created:  2022-07-06                                                                                                    #
# Last updated:  2022-07-06                                                                                                    #
# Source:        FishBT_2022                                                                                                   #
#                                                                                                                              #
# Comment:       Distribution of migration category, number of studies, sample size tested, and partial migration evidence     #
#                using biogeochemical tags in the phylogeny of global migratory fishes (Actinopterygii)                        #
#                                                                                                                              #
################################################################################################################################

############################### Packages loading ###############################
rm(list = ls())
options(warn = -1)

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(patchwork)
library(ggplot2)
library(cowplot)
library(ggimage)
library(tidyverse)
library(reshape2)
library(ggstar)
library(ggh4x)

################################# Data loading #################################
res_phylo <- readRDS("input/res_phylo.rds")
dt <- readRDS("input/tags_2022.03.21.rds")
phylopic_info <- readRDS("input/phylopic_info.rds")
fish <- readRDS("input/fish.rds")
partial_migration <- readRDS("input/partial_migration.rds")

################################ Data cleaning #################################
dt$N_analysis <- as.numeric(dt$N_analysis)

tag <- dt %>%
  dplyr::select(Binomial_fishbase, Tags_type, N_analysis) %>%
  group_by(Binomial_fishbase, Tags_type) %>%
  mutate(n = sum(N_analysis, na.rm = T)) %>%
  dplyr::select(-N_analysis) %>%
  distinct()
names(dt)

phylopic_info$groms <- factor(phylopic_info$groms,
  levels = c(
    "Amphidromous",
    "Anadromous",
    "Catadromous",
    "Oceanodromous",
    "Potamodromous"
  )
)

study <- dt %>%
  dplyr::select(Pdf_id, Binomial_fishbase, Element_type) %>%
  group_by(Binomial_fishbase, Element_type) %>%
  count()

fish$species$Species <- gsub(" ", "_", fish$species$Species)
data_fish <- res_phylo$Insertions_data %>%
  left_join(fish$species[, c(2, 26)], by = c("s" = "Species")) %>%
  dplyr::select(-insertions)
names(data_fish)[4] <- "GROMs"

data_fish$GROMs <- str_to_title(data_fish$GROMs)
data_fish$GROMs <- factor(data_fish$GROMs,
  levels = c(
    "Amphidromous",
    "Anadromous",
    "Catadromous",
    "Oceanodromous",
    "Potamodromous"
  )
)

################################### Ploting ####################################
p <- ggtree::ggtree(res_phylo$Phylogeny,
  layout = "fan",
  open.angle = 5,
  col = "grey70"
)
# GROMs_category
p <- p + new_scale_fill() +
  geom_fruit(
    data = data_fish,
    geom = geom_tile,
    mapping = aes(y = s, fill = GROMs), size = 0,
    width = 16,
    offset = 0.035,
    pwidth = 0.5
  ) +
  scale_fill_manual(
    name = "GROMS category", values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500"),
    guide = guide_legend(order = 3, reverse = T)
  )

# Number of studies
p <- p + new_scale_fill() +
  geom_fruit(
    data = study, geom = geom_bar,
    mapping = aes(x = n, y = Binomial_fishbase, fill = Tags_type),
    orientation = "y", offset = 0.035, pwidth = 0.15,
    stat = "identity",
    axis.params = list(
      title = "",
      axis = "x",
      title.size = 1.5,
      title.height = 0.01,
      text.size = 1.5,
      text.angle = 0,
      vjust = 1.3,
      nbreak = 3
    ), grid.params = list(linetype = 2)
  ) +
  scale_fill_manual(
    name = "Biogeochemical tag",
    values = c("#66C2A5", "#FC8D62"),
    guide = guide_legend(keywidth = 0.65, keyheight = 0.65, order = 1, reverse = T)
  )

# number of tagged
tag$logn <- log10(tag$n)
p <- p + new_scale_fill() +
  geom_fruit(
    data = tag, geom = geom_bar,
    mapping = aes(x = logn, y = Binomial_fishbase, fill = Tags_type),
    orientation = "y", offset = 0.04, pwidth = 0.15,
    stat = "identity",
    axis.params = list(
      title = "",
      axis = "x",
      title.size = 1.5,
      title.height = 0.01,
      text.size = 1.5,
      text.angle = 0,
      vjust = 1.3,
      nbreak = 3
    ), grid.params = list(linetype = 2)
  ) +
  scale_fill_manual(
    name = "Biogeochemical tag",
    values = c("#66C2A5", "#FC8D62"),
    guide = "none"
  )

# partial_migration
p <- p + new_scale_fill() +
  geom_fruit(
    data = partial_migration[which(partial_migration$Element_type == "Microchemistry"), ],
    geom = geom_point,
    mapping = aes(y = Binomial_fishbase, fill = Element_type),
    shape = 1,
    size = 1,
    position = position_identityx(),
    pwidth = 0.1,
    colour = "#FC8D62",
    offset = 0.03,
    show.legend = F
  )

p <- p + new_scale_fill() +
  geom_fruit(
    data = partial_migration[which(partial_migration$Element_type == "Isotope"), ],
    geom = geom_point,
    mapping = aes(y = Binomial_fishbase, fill = Element_type),
    shape = 1,
    size = 1,
    position = position_identityx(),
    pwidth = 0.1,
    colour = "#66C2A5",
    offset = 0.03,
    show.legend = F
  )

# phylopic_info
p <- p + new_scale_color() +
  geom_fruit(
    data = phylopic_info,
    geom = geom_phylopic,
    mapping = aes(y = taxon, image = uid, color = groms),
    size = 0.06,
    # color = "grey50",
    offset = 0.15,
    alpha = 0.8,
    position = position_identityx()
  ) +
  scale_color_manual(
    name = "GROMS category", values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500"),
    guide = guide_legend(order = 3, reverse = T)
  ) +
  theme(
    legend.title = element_text(size = 5),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing = unit(0.1, "cm"),
    legend.text = element_text(size = 5),
    legend.box.spacing = unit(-1, "cm")
  )

p <- p +
  annotate("text", x = 350, y = 0, label = "I", size = 1.5, vjust = 3) +
  annotate("text", x = 388, y = 0, label = "II", size = 1.5, vjust = 3) +
  annotate("text", x = 450, y = 0, label = "III", size = 1.5, vjust = 3) +
  annotate("text", x = 485, y = 0, label = "IV", size = 1.5, vjust = 3)

min_map <- ggplot(tag) +
  geom_density(aes(n, fill = Tags_type), col = NA, alpha = 0.75, adjust = 1) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62"), guide = "none") +
  labs(
    x = "Sample size tested",
    y = "Density"
  ) +
  scale_x_log10(
    guide = "axis_logticks",
    limits = c(.15, 30000),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c(
      expression(10^0),
      expression(10^1),
      expression(10^2),
      expression(10^3),
      expression(10^4)
    ),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.5), breaks = c(0, 0.25, 0.5)) +
  theme_classic() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text = element_text(size = 4, colour = "black"),
    axis.title = element_text(size = 5),
    axis.line = element_line(size = 0.1, colour = "black"),
    axis.ticks = element_line(size = 0.1, colour = "black")
  )

################################ Multiple plots ################################
ggdraw() +
  draw_plot(p, 0, 0, 1, 1) +
  draw_plot(min_map, 0.36, 0.42, 0.13, 0.13) +
  theme(plot.margin = margin(0, 0, 0, 0))

################################# Fig.3 saving #################################
ggsave(filename = "output/figure3.pdf", width = 18, height = 18, dpi = 600, units = "cm", family = "serif")
