################################################################################################################################
#                                                                                                                              #
# Purpose:       Fig.4 Plotting                                                                                                #
#                                                                                                                              #
# Author:        Liuyong Ding et al.                                                                                           #
# Contact:       ly_ding@126.com                                                                                               #
# Client:        Liuyong Ding et al.                                                                                           #
#                                                                                                                              #
# Code created:  2022-07-06                                                                                                    #
# Last updated:  2022-07-06                                                                                                    #
# Source:        FishBT_2022                                                                                                   #
#                                                                                                                              #
# Comment:       Proportion distribution of microchemical and isotopic elements tested using multiple tissues of migratory     #
#                fish across the periodic table of the elements                                                                #
#                                                                                                                              #
################################################################################################################################

############################### Packages loading ###############################
rm(list = ls())
options(warn = -1)

library(readxl)
library(ggplot2)
library(dplyr)
library(magrittr) # extract2()
library(ggnewscale)
library(ggimage)
library(magick)
library(grid)
library(ggstar)

################################# Data loading #################################
values <- readRDS("input/periodic_data.rds")
data <- readRDS("input/tags_2022.03.21.rds")
fish <- image_read("input/fish.svg")

################################ Data cleaning #################################
oto_m <- data %>%
  dplyr::filter(Tags_type == "Microchemistry" & Element_tissues == "Otolith") %>%
  dplyr::select(Element_composition)

oto_m <- unlist(strsplit(oto_m$Element_composition, split = ";")) %>%
  as.data.frame() %>%
  distinct() %>%
  set_names("oto_m")

oto_m <- unique(oto_m$oto_m)
oto_m <- oto_m[-which(nchar(oto_m) < 1)]
length(oto_m)

oto_i <- data %>%
  dplyr::filter(Tags_type == "Isotope" & Element_tissues == "Otolith") %>%
  dplyr::select(Element_composition)

oto_i <- unlist(strsplit(oto_i$Element_composition, split = ";")) %>%
  as.data.frame() %>%
  distinct() %>%
  set_names("oto_i")
oto_i$oto_i <- gsub("[^A-Z|^a-z]", "", oto_i$oto_i)
oto_i$oto_i <- ifelse(oto_i$oto_i == "dC", "C", oto_i$oto_i)
oto_i$oto_i <- ifelse(oto_i$oto_i == "dN", "N", oto_i$oto_i)
oto_i$oto_i <- ifelse(oto_i$oto_i == "dNd", "Nd", oto_i$oto_i)
oto_i$oto_i <- ifelse(oto_i$oto_i == "dO", "O", oto_i$oto_i)
oto_i$oto_i <- ifelse(oto_i$oto_i == "dS", "S", oto_i$oto_i)

oto_i <- unique(oto_i$oto_i)
length(oto_i)

oto_all <- unique(c(oto_m, oto_i))

eye <- data %>%
  dplyr::filter(Tags_type == "Isotope" & Element_tissues == "Eye_lense") %>%
  dplyr::select(Element_composition)

eye <- c("N", "C", "S")


data$Element_tissues <- ifelse(data$Element_tissues == "Vertebra", "Vertebrae", data$Element_tissues)

table(data$Element_tissues)

Vertebrae_i <- data %>%
  dplyr::filter(Tags_type == "Isotope" & Element_tissues == "Vertebrae") %>%
  dplyr::select(Element_composition) %>%
  distinct()

Vertebrae_i <- c("C", "N", "O", "S", "Sr", "Nd")

Vertebrae_m <- data %>%
  dplyr::filter(Tags_type == "Microchemistry" & Element_tissues == "Vertebrae") %>%
  dplyr::select(Element_composition)

Vertebrae_m <- unlist(strsplit(Vertebrae_m$Element_composition, split = ";")) %>%
  as.data.frame() %>%
  distinct() %>%
  set_names("Vertebrae_m")

Vertebrae_m <- unique(Vertebrae_m$Vertebrae_m)

Vertebrae_all <- unique(c(Vertebrae_m, Vertebrae_i))

Scale_m <- data %>%
  dplyr::filter(Tags_type == "Microchemistry" & Element_tissues == "Scale") %>%
  dplyr::select(Element_composition)
Scale_m <- unlist(strsplit(Scale_m$Element_composition, split = ";")) %>%
  as.data.frame() %>%
  distinct() %>%
  set_names("Scale")

Scale_m <- unique(Scale_m$Scale)

Scale_i <- data %>%
  dplyr::filter(Tags_type == "Isotope" & Element_tissues == "Scale") %>%
  dplyr::select(Element_composition) %>%
  distinct()
Scale_i <- c("C", "N", "H", "Sr", "S")

Scale_all <- unique(c(Scale_m, Scale_i))

dt1 <- subset(data, data$Tags_type == "Microchemistry")
# res <- dt1[5:8,]
res1 <- unlist(strsplit(dt1$Element_composition, split = ";"))
n.m <- data.frame(s1 = res1) %>%
  group_by(s1) %>%
  count() %>%
  set_names(c("Microchemistry", "n_m"))
n.m <- n.m[-1, ]

values$micro <- ifelse(values$Symbol %in% res1, 1, 0)
values <- left_join(values, n.m, by = c("Symbol" = "Microchemistry"))
values$n_m <- ifelse(values$micro == 1, values$n_m / sum(values$n_m, na.rm = T) * 100, NA)
values$n_m <- cut(values$n_m, breaks = seq(0, max(values$n_m, na.rm = T) + 4, by = 3))


dt2 <- subset(data, data$Tags_type == "Isotope")
res2 <- unlist(strsplit(dt2$Element_composition, split = ";"))

n.i <- data.frame(s2 = gsub("[^A-Z|^a-z]", "", res2)) %>%
  group_by(s2) %>%
  count() %>%
  set_names(c("Isotope", "n_i"))
n.i

n.i$Isotope <- ifelse(n.i$Isotope == "dC", "C", n.i$Isotope)
n.i$Isotope <- ifelse(n.i$Isotope == "dN", "N", n.i$Isotope)
n.i$Isotope <- ifelse(n.i$Isotope == "dNd", "Nd", n.i$Isotope)
n.i$Isotope <- ifelse(n.i$Isotope == "dO", "O", n.i$Isotope)
n.i$Isotope <- ifelse(n.i$Isotope == "dS", "S", n.i$Isotope)

n.i
values$iso <- ifelse(values$Symbol %in% n.i$Isotope, 1, 0)
values <- left_join(values, n.i, by = c("Symbol" = "Isotope"))
values$n_i <- ifelse(values$iso == 1, values$n_i / sum(values$n_i, na.rm = T) * 100, NA)

################################### Ploting ####################################
p <- ggplot() +
  # micro
  geom_point(
    data = values %>% filter(micro == 1),
    size = 16,
    shape = 15,
    # colour = "#a58394",
    aes(y = IUPAC_Period, x = IUPAC_Group, colour = n_m)
  ) +
  # iso
  # new_scale_colour()+
  geom_point(
    data = values %>% filter(iso > 0),
    # size = 1.5,
    shape = 1,
    colour = "black",
    aes(y = IUPAC_Period + 0.3, x = IUPAC_Group, size = n_i)
  ) +
  geom_star(
    data = values[which(values$Symbol %in% oto_m), ] %>%
      dplyr::select(IUPAC_Period, IUPAC_Group, Symbol) %>%
      distinct(),
    aes(y = IUPAC_Period - 0.36, x = IUPAC_Group + 0.36), starshape = 11,
    size = 1.6, fill = "grey50", colour = NA
  ) + 
  geom_star(
    data = values[which(values$Symbol %in% eye), ] %>%
      dplyr::select(IUPAC_Period, IUPAC_Group, Symbol) %>%
      distinct(),
    aes(y = IUPAC_Period - 0.12, x = IUPAC_Group + 0.36), starshape = 15,
    size = 1.6, fill = "grey50", colour = NA
  ) + 
  geom_star(
    data = values[which(values$Symbol %in% Vertebrae_all), ] %>%
      dplyr::select(IUPAC_Period, IUPAC_Group, Symbol) %>%
      distinct(),
    aes(y = IUPAC_Period + 0.12, x = IUPAC_Group + 0.36), starshape = 22,
    size = 1.6, fill = "grey50", colour = NA
  ) + 
  geom_star(
    data = values[which(values$Symbol %in% Scale_all), ] %>%
      dplyr::select(IUPAC_Period, IUPAC_Group, Symbol) %>%
      distinct(),
    aes(y = IUPAC_Period + 0.36, x = IUPAC_Group + 0.36), starshape = 4,
    size = 1.6, fill = "grey50", colour = NA
  ) + 
  # boxes for all elements
  geom_point(
    data = values,
    size = 16,
    shape = 0,
    aes(
      y = IUPAC_Period,
      x = IUPAC_Group
    )
  ) +
  # element symbol
  geom_text(
    data = values,
    size = 3,
    colour = "black",
    fontface = "bold",
    aes(
      label = Symbol,
      y = IUPAC_Period - 0.2, x = IUPAC_Group
    )
  ) +
  # atomic number
  geom_text(
    data = values,
    size = 1.8,
    colour = "black",
    aes(
      label = IUPAC_Number,
      y = IUPAC_Period - 0.32,
      x = ifelse(IUPAC_Number %in% c("57-71", "89-103"), IUPAC_Group, IUPAC_Group - 0.32)
    )
  ) +
  # element name
  geom_text(
    data = values,
    size = 1.6,
    colour = "black",
    aes(
      label = tolower(Name),
      y = IUPAC_Period + 0.05, x = IUPAC_Group
    )
  ) +
  geom_text(
    data = values %>% filter(IUPAC_Period < 8),
    size = 1.7,
    colour = "black",
    aes(
      label = IUPAC_Period,
      y = IUPAC_Period
    ),
    x = 0.38
  ) +
  # group labels # positions manually adjusted
  geom_text(
    data =
      data.frame(
        y = c(0.38, 1.38, rep(3.38, 10), rep(1.38, 5), 0.38),
        x = seq(1, 18)
      ),
    size = 1.7,
    colour = "black",
    aes(label = x, x = x, y = y)
  ) +
  scale_x_continuous(
    breaks = seq(
      min(values$IUPAC_Group),
      max(values$IUPAC_Group)
    ),
    limits = c(
      min(values$IUPAC_Group) - 1,
      max(values$IUPAC_Group) + 1
    ),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = "reverse",
    breaks = seq(
      min(values$IUPAC_Period),
      max(values$IUPAC_Period)
    ),
    limits = c(
      max(values$IUPAC_Period) + 1,
      min(values$IUPAC_Period) - 1.5
    ),
    expand = c(0, 0)
  ) +
  # coord_equal()+
  scale_size_continuous("Isotope (%)", range = c(2, 5)) +
  scale_colour_viridis_d("Microchemistry (%)",
    option = "inferno",
    direction = -1,
    na.value = "grey70",
    end = 0.9,
    begin = 0.2
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.25, 0.90),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(0.2, "line"),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  guides(
    colour = guide_legend(nrow = 2, override.aes = list(size = 3), title.position = "top", order = 1),
    size = guide_legend(nrow = 1, title.position = "top", order = 2)
  )

############################### Legends Plotting ###############################
lb <- ggplot(data = data.frame(x = 0, y = 0), aes(x, y)) +
  geom_point(size = 16, shape = 0) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  annotate("text", x = 0, y = 0.05, label = "Sr", size = 3, fontface = "bold") +
  annotate("text", x = 0, y = -0.05, label = "strontium", size = 2) +
  annotate("text", x = -0.10, y = 0.15, label = "38", size = 2) +
  annotate("segment", x = -0.3, xend = -0.16, y = 0.15, yend = 0.15, size = 0.4) +
  annotate("text", x = -0.57, y = 0.15, label = "Atomic number", size = 2) +
  annotate("segment", x = -0.3, xend = -0.08, y = 0.05, yend = 0.05, size = 0.4) +
  annotate("text", x = -0.45, y = 0.05, label = "Symbol", size = 2) +
  annotate("segment", x = -0.3, xend = -0.17, y = -0.05, yend = -0.05, size = 0.4) +
  annotate("text", x = -0.43, y = -0.05, label = "Name", size = 2) +
  theme_void() +
  coord_equal()


p +
  annotation_custom(ggplotGrob(lb), xmin = -7, xmax = 12, ymin = -7, ymax = -12) +
  annotation_custom(rasterGrob(image_resize(fish, 900)), xmin = 7, xmax = 12, ymin = -7, ymax = 2.5) +
  geom_star(
    data = data.frame(x = 8.95, y = 1.3),
    aes(x = x, y = y), starshape = 11,
    size = 3, fill = "grey50", colour = NA
  ) +
  geom_star(
    data = data.frame(x = 7.75, y = 1.35),
    aes(x = x, y = y), starshape = 15,
    size = 3, fill = "grey50", colour = NA
  ) +
  geom_star(
    data = data.frame(x = 8.5, y = 2.95),
    aes(x = x, y = y), starshape = 22,
    size = 3, fill = "grey50", colour = NA
  ) +
  geom_star(
    data = data.frame(x = 10.85, y = 1.35),
    aes(x = x, y = y), starshape = 4,
    size = 3, fill = "grey50", colour = NA
  )

################################# Fig.4 saving #################################
# size 16 leaves no gaps between the element boxes at fig.width=9.51,fig.height=5.51
ggsave(filename = "output/figure4.pdf", width = 9.51, height = 5.51, dpi = 600)
