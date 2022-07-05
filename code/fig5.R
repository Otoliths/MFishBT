################################################################################################################################
#                                                                                                                              #
# Purpose:       Fig.5 Plotting                                                                                                #
#                                                                                                                              #
# Author:        Liuyong Ding et al.                                                                                           #
# Contact:       ly_ding@126.com                                                                                               #
# Client:        Liuyong Ding et al.                                                                                           #
#                                                                                                                              #
# Code created:  2022-07-06                                                                                                    #
# Last updated:  2022-07-06                                                                                                    #
# Source:        FishBT_2022                                                                                                   #
#                                                                                                                              #
# Comment:       Sample size and value distribution of major microchemical elements and isotopic ratio detected in otolith     #
#                core areas of migratory fishes across GROMS category                                                          #
#                                                                                                                              #
################################################################################################################################

############################### Packages loading ###############################
rm(list = ls())
options(warn = -1)

library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)
library(readxl)
library(stringr)
library(SDMTools)
library(ggforce)
library(ggtext)
library(ggbreak)
library(ggh4x)
library(rcompanion)
library(nparcomp)

################################# Data loading #################################
oto <- readRDS("input/otolith_core.rds")


################################ Data cleaning #################################
d1 <- oto$oto.core.Sr.Ca

n1 <- d1 %>%
  dplyr::select(GROMs_fishbase, Mean_fixed, N) %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N)) %>%
  mutate(label = paste0("n=", n), x = rep(0.02, 5))
n1


point1 <- d1 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

# Kruskal-Wallis
kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d1)

# r1 <- npmc(Mean_fixed~GROMs_fishbase, d1,method = "holm")
# summary(r1)
# plot(r1)

d1.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = d1, asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d1.mc$connames, p.value = d1.mc$Analysis$p.Value)
d1.mc$Data.Info

n1 <- n1 %>% mutate(npm = c("d", "b", "e", "c", "a"))

################################################################################
#                                                                              #
#                           oto.core Sr/Ca plotting                            #
#                                                                              #
################################################################################

p1 <- ggplot(
  data = d1,
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, fill = NA, show.legend = F) +
  geom_text(data = n1, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n1, aes(y = x - 0.01, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 0, col = "grey60", fontface = "italic") +
  geom_point(
    data = point1,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point1,show.legend = F,lty = 2,size = 0.3)+
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 0.02, label = "***", colour = "grey70", size = 2, hjust = 0) +
  theme_bw() +
  coord_flip() +
  ylab("Sr/Ca") +
  xlab("") +
  scale_y_log10(
    guide = "axis_logticks",
    limits = c(0.0000001, 0.1),
    breaks = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1),
    labels = c(
      expression(10^-7),
      expression(10^-6),
      expression(10^-5),
      expression(10^-4),
      expression(10^-3),
      expression(10^-2),
      expression(10^-1)
    )
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6, colour = "black"),
    panel.border = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    ggh4x.axis.ticks.length.minor = rel(0.4),
    ggh4x.axis.ticks.length.mini = rel(0.3)
  )


################################ Data cleaning #################################
d2 <- oto$oto.core.Sr

s <- data.frame(GROMs_fishbase = "Oceanodromous", Mean_fixed = NA, N = 0)

n2 <- rbind(d2, s) %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N, na.rm = T)) %>%
  mutate(label = paste0("n=", n), x = rep(0.765, 5))
# n2[4,3] <- NA
n2

point2 <- d2 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

# Kruskal-Wallis
kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d2)

# r2 <- npmc(Mean_fixed~GROMs_fishbase, d2,method = "holm")
# summary(r2)
# plot(r2)

d2.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = d2, asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d2.mc$connames, p.value = d2.mc$Analysis$p.Value)
d2.mc$Data.Info

n2 <- n2 %>% mutate(npm = c("bc", "a", "b", NA, "c"))

################################################################################
#                                                                              #
#                             oto.core Sr plotting                             #
#                                                                              #
################################################################################

p2 <- ggplot(
  data = rbind(d2, s),
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, show.legend = F, fill = NA) +
  geom_text(data = n2, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n2, aes(y = x - 0.01, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 0, col = "grey60", fontface = "italic") +
  geom_point(
    data = point2,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point2,show.legend = F,lty = 2,size = 0.3)+
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 0.765, label = "***", colour = "grey70", size = 2, hjust = 0) +
  theme_bw() +
  coord_flip() +
  ylab(expression(""^87 * "Sr/"^86 * "Sr" * "(‰)")) +
  xlab("") +
  scale_y_continuous(
    expand = c(0.001, 0.001),
    limits = c(0.703, 0.775),
    breaks = c(0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77)
  ) +
  # scale_x_discrete(expand = c(0.001,0.003))+
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.border = element_blank()
  )



################################ Data cleaning #################################
d3 <- oto$oto.core.Mg.Ca

n3 <- d3 %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N, na.rm = T)) %>%
  mutate(label = paste0("n=", n), x = rep(0.02, 5))


point3 <- d3 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

point3[4, 4] <- 0

kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d3)

# r3 <- npmc(Mean_fixed~GROMs_fishbase, d3,method = "holm")
# summary(r3)
# plot(r3)


d3.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = d3, asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d3.mc$connames, p.value = d3.mc$Analysis$p.Value)
d3.mc$Data.Info


n3 <- n3 %>% mutate(npm = c("ab", "a", "ab", "b", "ab"))

################################################################################
#                                                                              #
#                           oto.core Mg/Ca plotting                            #
#                                                                              #
################################################################################

p3 <- ggplot(
  data = d3,
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, show.legend = F, fill = NA) +
  geom_point(
    data = point3,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point3,show.legend = F,lty = 2,size = 0.3)+
  geom_text(data = n3, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n3, aes(y = x - 0.01, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 1, col = "grey60", fontface = "italic") +
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 0.02, label = "***", colour = "grey70", size = 2, hjust = 0) +
  theme_bw() +
  coord_flip() +
  ylab("Mg/Ca") +
  xlab("") +
  scale_y_log10(
    guide = "axis_logticks",
    limits = c(0.0000001, 0.1),
    breaks = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1),
    labels = c(
      expression(10^-7),
      expression(10^-6),
      expression(10^-5),
      expression(10^-4),
      expression(10^-3),
      expression(10^-2),
      expression(10^-1)
    )
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.border = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    ggh4x.axis.ticks.length.minor = rel(0.4),
    ggh4x.axis.ticks.length.mini = rel(0.3)
  )

p3 <- p3 + theme(axis.text.y = element_blank())

################################ Data cleaning #################################
d4 <- oto$oto.core.O

n4 <- d4 %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N, na.rm = T)) %>%
  mutate(label = paste0("n=", n), x = rep(3, 5))


point4 <- d4 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

# Kruskal-Wallis
kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d4)

# r4 <- npmc(Mean_fixed~GROMs_fishbase, subset(d4,d4$GROMs_fishbase!= "Amphidromous"),method = "holm")
# summary(r4)
# plot(r4)

d4.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = subset(d4, d4$GROMs_fishbase != "Amphidromous"), asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d4.mc$connames, p.value = d4.mc$Analysis$p.Value)
d4.mc$Data.Info


n4 <- n4 %>% mutate(npm = c(NA, "c", "a", "b", "a"))

################################################################################
#                                                                              #
#                             oto.core O plotting                              #
#                                                                              #
################################################################################

p4 <- ggplot(
  data = d4,
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, show.legend = F, fill = NA) +
  geom_text(data = n4, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n4, aes(y = x - 2, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 1, col = "grey60", fontface = "italic") +
  geom_point(
    data = point4,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point4,show.legend = F,lty = 2,size = 0.3)+
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 3, label = "***", colour = "grey70", size = 2, hjust = 0) +
  theme_bw() +
  coord_flip() +
  ylab(expression(delta * ""^18 * "O" * "(‰)")) +
  xlab("") +
  scale_y_continuous(
    expand = c(0.001, 0.001),
    limits = c(-14.8, 6),
    breaks = c(-15, -10, -5, 0, 5)
  ) +
  # scale_x_discrete(expand = c(0.001,0.003))+
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.border = element_blank()
  )

p4 <- p4 + theme(axis.text.y = element_blank())


################################ Data cleaning #################################
d5 <- oto$oto.core.Ba.Ca

n5 <- d5 %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N, na.rm = T)) %>%
  mutate(label = paste0("n=", n), x = rep(0.02, 5))


point5 <- d5 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

point5[4, 4] <- 0
# Kruskal-Wallis
kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d4)

# r5 <- npmc(Mean_fixed~GROMs_fishbase, d5,method = "holm")
# summary(r5)
# plot(r5)

d5.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = d5, asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d5.mc$connames, p.value = d5.mc$Analysis$p.Value)
d5.mc$Data.Info

n5 <- n5 %>% mutate(npm = c("b", "a", "b", "c", "d"))

################################################################################
#                                                                              #
#                           oto.core Ba/Ca plotting                            #
#                                                                              #
################################################################################

p5 <- ggplot(
  data = d5,
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, show.legend = F, fill = NA) +
  geom_point(
    data = point5,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point5,show.legend = F,lty = 2,size = 0.3)+
  geom_text(data = n5, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n5, aes(y = x - 0.01, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 1, col = "grey60", fontface = "italic") +
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 0.02, label = "***", colour = "grey70", size = 2, hjust = 0) +
  theme_bw() +
  coord_flip() +
  ylab("Ba/Ca") +
  xlab("") +
  scale_y_log10(
    guide = "axis_logticks",
    limits = c(0.0000001, 0.1),
    breaks = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1),
    labels = c(
      expression(10^-7),
      expression(10^-6),
      expression(10^-5),
      expression(10^-4),
      expression(10^-3),
      expression(10^-2),
      expression(10^-1)
    )
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.border = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    ggh4x.axis.ticks.length.minor = rel(0.4),
    ggh4x.axis.ticks.length.mini = rel(0.3)
  )
p5 <- p5 + theme(axis.text.y = element_blank())

################################ Data cleaning #################################
d6 <- oto$oto.core.C

ss <- data.frame(GROMs_fishbase = "Catadromous", Mean_fixed = NA, N = 0)
d6 <- rbind(d6, ss)
n6 <- d6 %>%
  group_by(GROMs_fishbase) %>%
  summarise(n = sum(N, na.rm = T)) %>%
  mutate(label = paste0("n=", n), x = rep(1.4, 5))


point6 <- d6 %>%
  group_by(GROMs_fishbase) %>%
  summarise(
    wt_sd = wt.sd(Mean_fixed, N),
    wt_mean = wt.mean(Mean_fixed, N)
  ) %>%
  mutate(
    lower = wt_mean - wt_sd,
    upper = wt_mean + wt_sd,
    Mean_fixed = wt_mean
  )

# Kruskal-Wallis
kruskal.test(Mean_fixed ~ GROMs_fishbase, data = d6)

# r6 <- npmc(Mean_fixed~GROMs_fishbase, d6[which(d6$GROMs_fishbase %in% c("Potamodromous","Oceanodromous","Anadromous")),],method = "holm")
# summary(r6)
# plot(r6)

d6.mc <- mctp(Mean_fixed ~ GROMs_fishbase,
  data = subset(d6, d6$GROMs_fishbase != "Amphidromous"), asy.method = "mult.t",
  type = "Tukey", alternative = "two.sided",
  plot.simci = F, info = FALSE
)
cldList(comparison = d6.mc$connames, p.value = d6.mc$Analysis$p.Value)
d6.mc$Data.Info

n6 <- n6 %>% mutate(npm = c(NA, "b", NA, "a", "ab"))

################################################################################
#                                                                              #
#                             oto.core C plotting                              #
#                                                                              #
################################################################################

p6 <- ggplot(
  data = d6,
  aes(y = Mean_fixed, x = GROMs_fishbase)
) +
  geom_sina(aes(size = N, col = GROMs_fishbase), alpha = 0.1, fill = NA, show.legend = F) +
  geom_text(data = n6, aes(y = x, x = GROMs_fishbase, label = label, colour = GROMs_fishbase), show.legend = F, size = 2, hjust = 0) +
  geom_text(data = n6, aes(y = x - 3, x = GROMs_fishbase, label = npm), show.legend = F, size = 2, hjust = 0, col = "grey60", fontface = "italic") +
  geom_point(
    data = point6,
    aes(y = wt_mean, x = GROMs_fishbase), show.legend = F, col = "grey60", size = 1
  ) +
  # geom_segment(aes(y = lower,x = GROMs_fishbase,yend = upper, xend = GROMs_fishbase),col = "grey60", data = point6,show.legend = F,lty = 2,size = 0.3)+
  scale_fill_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  scale_colour_manual(values = c("#b30059", "#220050", "#0091a8", "#359023", "#ffa500")) +
  annotate("text", x = 5.2, y = 1.4, label = "***", colour = "grey70", size = 2, hjust = 0) +
  # scale_size_continuous("Size",range = c(0.5,2))+
  theme_bw() +
  coord_flip() +
  ylab(expression(delta * ""^13 * "C" * "(‰)")) +
  xlab("") +
  scale_y_continuous(
    expand = c(0.001, 0.001),
    limits = c(-26, 6),
    breaks = c(-25, -20, -15, -10, -5, 0, 5)
  ) +
  # scale_x_discrete(expand = c(0.001,0.003))+
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.border = element_blank()
  )

p6 <- p6 + theme(axis.text.y = element_blank())


################################################################################
#                                                                              #
# Silhouette images of fishes plotting,these images of fishes were sourced     #
# from http://www.phylopic.org/.                                               #
#                                                                              #
################################################################################
phylopic_info_m <- data.frame(
  GROMs_fishbase = c(
    "Amphidromous", # Acanthogobius_flavimanus
    "Anadromous", # Mugil_cephalus
    "Catadromous", # Anguilla_anguilla
    "Oceanodromous", # Thunnus_albacares
    "Potamodromous"
  ), # Cyprinus_carpio
  uid = c(
    "966455e3-87b9-4bb3-86ef-7856e8cc3456",
    "ead318d0-2ded-4ca5-b263-b9702374a564",
    "629e1135-d0c5-4f02-bb13-8cbf63e48692",
    "be13bcdf-ebef-4288-bdc6-114751cdb550",
    "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1"
  )
)


labels_m <- setNames(
  paste0(
    "<img src='http://www.phylopic.org/assets/images/submissions/",
    phylopic_info_m$uid, ".512.png'  width='20' /><br>",
    sapply(
      strwrap(phylopic_info_m$GROMs_fishbase, width = 10, simplify = FALSE),
      function(x) paste(x, collapse = "<br>")
    )
  ),
  phylopic_info_m$GROMs_fishbase
)


labels_m[5] <- "<img src='input/phylopic/25a76141-cca2-4b56-bd9f-9f5fdcd8efa1.512.png'  width='20' /><br>Potamodromous"


phylopic_info_i <- data.frame(
  GROMs_fishbase = c(
    "Amphidromous", # Gymnocypris_przewalskii
    "Anadromous", # Oncorhynchus_tshawytscha
    "Catadromous", # Lates_calcarifer
    "Oceanodromous", # Macruronus_magellanicus
    "Potamodromous"
  ), # Macquaria_ambigua
  uid = c(
    "ead318d0-2ded-4ca5-b263-b9702374a564",
    "9150be88-6910-4374-aa54-a7f8f3d79fb6",
    "08915a93-1ea3-4373-8ff0-72d05324537f",
    "8d92b454-3131-4bbd-ac9c-e1df14c2fc5a",
    "0dc8cdc8-eb45-4d23-aef4-9314968a23fb"
  )
)


labels_i <- setNames(
  paste0(
    "<img src='http://www.phylopic.org/assets/images/submissions/",
    phylopic_info_i$uid, ".512.png'  width='20' /><br>",
    sapply(
      strwrap(phylopic_info_i$GROMs_fishbase, width = 10, simplify = FALSE),
      function(x) paste(x, collapse = "<br>")
    )
  ),
  phylopic_info_i$GROMs_fishbase
)


p1 <- p1 + scale_x_discrete(labels = labels_m) +
  theme(axis.text.y = element_markdown(size = 6, colour = "black"))

p2 <- p2 + scale_x_discrete(labels = labels_i) +
  theme(axis.text.y = element_markdown(size = 6, colour = "black"))

################################ Multiple plots ################################
p1 + p3 + p5 + p2 + p6 + p4 + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "a") &
  theme(
    axis.text.x = element_text(size = 6, colour = "black"),
    panel.spacing = unit(0, "cm")
  )

################################# Fig.5 saving #################################
ggsave(filename = "output/figure5.pdf", width = 20, height = 12, dpi = 600, units = "cm")
