################################################################################################################################
#                                                                                                                              #
# Purpose:       Fig.1 Plotting                                                                                                #
#                                                                                                                              #
# Author:        Liuyong Ding et al.                                                                                           #
# Contact:       ly_ding@126.com                                                                                               #
# Client:        Liuyong Ding et al.                                                                                           #
#                                                                                                                              #
# Code created:  2022-07-06                                                                                                    #
# Last updated:  2022-07-06                                                                                                    #
# Source:        FishBT_2022                                                                                                   #
#                                                                                                                              #
# Comment:       Spatiotemporal distribution of migratory fish studies using biogeochemical tags                               #
#                                                                                                                              #
################################################################################################################################

############################### Packages loading ###############################
rm(list = ls())
options(warn = -1)
library(sf)
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(tidyverse)
library(ggspatial)
library(rworldxtra)
library(RColorBrewer)
library(classInt)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
library(ggrastr)
library(cowplot)
library(biscale)
library(geomtextpath)
library(dplyr)
library(tmap)
source("code/zzz.r")

################################# Data loading #################################
data("rivers")
fao <- readRDS("input/fao.rds")
point <- readRDS("input/occurence_point_2022.03.21.rds")
tags <- readRDS("input/tags_2022.03.21.rds")

################################ Data cleaning #################################
dt <- left_join(data_frame(
  Tags_type = rep(names(table(tags$Tags_type)),
    each = length(min(tags$Publication_date, na.rm = TRUE):max(tags$Publication_date, na.rm = TRUE))
  ),
  Publication_date = rep(
    min(tags$Publication_date, na.rm = TRUE):max(tags$Publication_date, na.rm = TRUE),
    length(names(table(tags$Tags_type)))
  )
),
tags[, c(2, 3, 16)] %>%
  distinct() %>%
  group_by(Tags_type) %>%
  count(Publication_date),
by = c("Tags_type", "Publication_date")
)
dt$n <- ifelse(is.na(dtt$n), 0, dtt$n)

############################### Fig1.b1 Plotting ###############################
b1 <- dt %>%
  group_by(Tags_type) %>%
  mutate(n = cumsum(n)) %>%
  ggplot(aes(x = Publication_date, y = n)) +
  geom_line(aes(col = Tags_type), size = 0.5, na.rm = T) +
  annotate("text", x = 2022, y = 666, label = "Microchemistry", size = 1.8, hjust = -0.1, vjust = -0.2, col = "#FC8D62") +
  annotate("text", x = 2022, y = 358, label = "Isotope", size = 1.8, hjust = -0.1, vjust = -0.2, col = "#66C2A5") +
  scale_color_manual(values = c("#66C2A5", "#FC8D62")) +
  labs(x = "Year", y = "Total no. of studies") + # Cumulative number of studies
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(1970, 2030), breaks = c(1970, 1980, 1990, 2000, 2010, 2020)) +
  scale_y_continuous(expand = c(0.02, 0.02), limits = c(0, 700)) +
  coord_fixed(0.04) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    # panel.border = element_blank(),
    legend.position = "none",
    # plot.margin = margin(0,3,0,0,unit = "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_line(size = 0.3, colour = "black"),
    axis.line = element_line(linetype = 1, color = "black", size = 0.3),
    axis.text.x = element_text(vjust = 0, colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title = element_text(size = 5, face = "bold"),
    plot.title = element_text(size = 6, hjust = 0, vjust = -3, face = "bold")
  )
b1

############################### Fig1.b2 Plotting ###############################
b2 <- point %>%
  group_by(Ecosystem) %>%
  count() %>%
  ggplot(aes(x = Ecosystem, y = n, fill = Ecosystem)) +
  geom_col(position = "dodge") +
  scale_fill_manual("Sampling point", values = c("#FFD700", "#1E90FF")) +
  labs(x = "", y = "No. of points") + # Cumulative number of studies
  theme_classic() +
  scale_x_discrete(expand = c(0.7, 0.05), label = c("Inland", "Marine")) +
  scale_y_continuous(expand = c(0.02, 0.02), limits = c(0, 6000)) +
  # coord_equal(0.0002)+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    # panel.border = element_blank(),
    legend.position = "none",
    # plot.margin = margin(0,3,0,0,unit = "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_line(size = 0.3, colour = "black"),
    axis.line = element_line(linetype = 1, color = "black", size = 0.3),
    axis.text.x = element_text(vjust = 0, colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title = element_text(size = 5, face = "bold"),
    plot.title = element_text(size = 6, hjust = 0, vjust = -3, face = "bold")
  )
b2

################################ Data cleaning #################################
iso <- point %>%
  filter(Element_type == "Isotope") %>%
  dplyr::select(Pdf_id, FAO_code, FAO_zone, Ecosystem) %>%
  distinct() %>%
  count(FAO_code) %>%
  rename(Isotope = n)

micro <- point %>%
  filter(Element_type == "Microchemistry") %>%
  dplyr::select(Pdf_id, FAO_code, FAO_zone, Ecosystem) %>%
  distinct() %>%
  count(FAO_code) %>%
  rename(Microchemistry = n)

fao_map <- fao %>%
  left_join(iso, by = c("code" = "FAO_code")) %>%
  left_join(micro, by = c("code" = "FAO_code"))

pnts_sf <- st_as_sf(point, coords = c("Longitude", "Latitude"), crs = st_crs(fao_map))
dt1 <- bi_class2(fao_map[1:7, ], x = Microchemistry, y = Isotope, style = "quantile", dim = 4)
dt2 <- bi_class2(fao_map[8:26, ], x = Microchemistry, y = Isotope, style = "quantile", dim = 4)

############################### Legends Plotting ###############################
# Now we need to create a color matrix. Copy and paste the following function into R:
colmat <- function(nquantiles = 10,
                   upperleft = rgb(0, 150, 235, maxColorValue = 255),
                   upperright = rgb(130, 0, 80, maxColorValue = 255),
                   bottomleft = "grey",
                   bottomright = rgb(255, 230, 15, maxColorValue = 255),
                   xlab = "x label", ylab = "y label") {
  my.data <- seq(0, 1, .01)
  my.class <- classIntervals(my.data, n = nquantiles, style = "quantile")
  my.pal.1 <- findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)
  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i, ] <- findColours(my.class, my.col)
  }
  plot(c(1, 1), pch = 19, col = my.pal.1, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1), frame.plot = F, xlab = xlab, ylab = ylab, cex.lab = 1.3)
  for (i in 1:101) {
    col.temp <- col.matrix[i - 1, ]
    points(my.data, rep((i - 1) / 100, 101), pch = 15, col = col.temp, cex = 1)
  }
  seqs <- seq(0, 100, (100 / nquantiles))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
}


# You can specify the number of quantiles, colors and labels of your color matrix. Example:
x1 <- colmat(
  nquantiles = 4,
  upperleft = "#73AE80",
  upperright = "#2A5A5B",
  bottomleft = "#ADD8E6",
  bottomright = "#6C83B5"
)
x1
m <- c(
  "4-4" = "#2A5A5B", "3-4" = "#427667", "2-4" = "#5A9273", "1-4" = "#73AE80",
  "4-3" = "#406779", "3-3" = "#578386", "2-3" = "#6E9F94", "1-3" = "#86BCA2",
  "4-2" = "#557597", "3-2" = "#6B91A5", "2-2" = "#82ADB4", "1-2" = "#99CAC3",
  "4-1" = "#6C83B5", "3-1" = "#819FC5", "2-1" = "#97BBD5", "1-1" = "#ADD8E6"
)
bi_legend2(m, dim = 4)


quantile(fao_map[8:26, ]$Isotope, na.rm = T)
quantile(fao_map[8:26, ]$Microchemistry, na.rm = T)

Marine_areas <- bi_legend2(m,
  dim = 4, xlab = "Microchemistry",
  ylab = "Isotope",
  size = 5
) +
  xlim(c("6", "19", "27", "101")) +
  ylim(c("4", "10", "21", "48")) +
  # scale_x_discrete(expand = c(0,0),labels = letters[1:4])+
  # scale_y_discrete(expand = c(0,0),labels = letters[1:4])+
  theme_bw() +
  labs(title = "Marine") +
  theme(
    panel.background = element_rect(fill = NA),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(vjust = 3, colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title = element_text(size = 5, face = "bold"),
    plot.title = element_text(size = 5, hjust = 0.5, vjust = -3, face = "bold")
  )
Marine_areas

x2 <- colmat(
  nquantiles = 4,
  upperleft = "#9972AF",
  upperright = "#AE3A4E",
  bottomleft = "#CABED0",
  bottomright = "#AF8E53"
)
x2 <- x2[-1, -1]
f <- c(
  "4-4" = x2[4, 4], "3-4" = x2[3, 4], "2-4" = x2[2, 4], "1-4" = x2[1, 4],
  "4-3" = x2[4, 3], "3-3" = x2[3, 3], "2-3" = x2[2, 3], "1-3" = x2[1, 3],
  "4-2" = x2[4, 2], "3-2" = x2[3, 2], "2-2" = x2[2, 2], "1-2" = x2[1, 2],
  "4-1" = x2[4, 1], "3-1" = x2[3, 1], "2-1" = x2[2, 1], "1-1" = x2[1, 1]
)

quantile(fao_map[1:7, ]$Isotope, na.rm = T)
quantile(fao_map[1:7, ]$Microchemistry, na.rm = T)

Inland_waters <- bi_legend2(f,
  dim = 4, xlab = "Microchemistry",
  ylab = "Isotope",
  size = 5
) +
  xlim(c("38", "91", "117", "178")) +
  ylim(c("17", "37", "54", "115")) +
  # scale_x_discrete(expand = c(0,0),labels = letters[1:4])+
  # scale_y_discrete(expand = c(0,0),labels = letters[1:4])+
  theme_bw() +
  labs(title = "Inland") +
  theme(
    panel.background = element_rect(fill = NA),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(vjust = 3, colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title = element_text(size = 5, face = "bold"),
    plot.title = element_text(size = 5, hjust = 0.5, vjust = -3, face = "bold")
  )
Inland_waters

################################ Fig1.a Plotting ###############################
ggplot() +
  geom_sf(data = dt1, mapping = aes(fill = bi_class), color = "white", size = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = f, na.value = "lightgrey") +
  new_scale_fill() +
  geom_sf(data = dt2, mapping = aes(fill = bi_class), color = "white", size = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = m, na.value = "lightgrey") +
  new_scale_colour() +
  geom_sf(data = rivers, aes(size = scalerank), color = "lightblue", show.legend = FALSE) +
  scale_size(range = c(0.01, 0.15)) +
  new_scale_colour() +
  geom_sf(data = pnts_sf, aes(colour = Ecosystem), size = 0.5, alpha = 1 / 2, fill = "black") +
  geom_sf_text(data = fao_map, aes(label = code), size = 2.8) +
  scale_color_manual("Sampling point", values = c("#FFD700", "#1E90FF")) + # c("#C8B35A","#4885C1")
  annotate("text", x = 182, y = 84, label = expression(bold("FAO major fishing areas")), size = 2, hjust = 0, vjust = -0.2) +
  annotate("text", x = 182, y = seq(-100, 80, length = 26), label = sort(paste(fao_map$code, fao_map$zone), decreasing = T), size = 2, hjust = 0) +
  annotate("text", x = -175, y = 85, label = "a", size = 3, fontface = "bold") +
  annotate("text", x = -175, y = -94, label = "b", size = 3, fontface = "bold") +
  annotate("text", x = -30, y = -96, label = "Legends", size = 2, fontface = "bold") +
  coord_sf(crs = 4326) +
  scale_x_continuous(limits = c(-180, 260), expand = c(0.01, 0.01)) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0.05, 0, 0.05, unit = "cm"),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.title = element_text(face = "bold")
  ) +
  annotation_custom(grob = ggplotGrob(Inland_waters), xmin = 55, xmax = 110, ymin = -230, ymax = 0) +
  annotation_custom(grob = ggplotGrob(Marine_areas), xmin = 120, xmax = 172.5, ymin = -230, ymax = 0) +
  annotation_custom(grob = ggplotGrob(b2), xmin = -15, xmax = 45, ymin = -145, ymax = -93) +
  annotation_custom(grob = ggplotGrob(b1), xmin = -180, xmax = -70, ymin = -145, ymax = -93)

################################# Fig.1 saving #################################
ggsave(filename = "output/figure1.pdf", width = 20, height = 14, dpi = 600, units = "cm")
