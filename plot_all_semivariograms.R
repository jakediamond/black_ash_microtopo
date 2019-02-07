# 
# Author: Jake Diamond
# Purpose: Plot all site semivariograms based on best detrends
# Date: February 5, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(ggpubr)

# Plot all data based on best detrend
# L3 no detrend
# L2 linear
# L1 no detrend
# D1 quad (linear??)
# D2 linear
# D3 no detrend
# D4 quad
# T1 quad
# T2 no detrend
# T3 no detrend
nd <- c("L3", "L1", "D3", "T2", "T3")
l <- c("L2", "D2")
q <- c("D1", "D4", "T1")

# Get all data
df <- read.csv("semi_results_quad.csv",
               stringsAsFactors = FALSE) %>%
  dplyr::filter(site %in% q) %>%
  bind_rows(read.csv("semi_results_linear.csv",
                     stringsAsFactors = FALSE) %>%
              dplyr::filter(site %in% l)) %>%
  bind_rows(read.csv("semi_results_no_trend.csv",
                     stringsAsFactor = FALSE) %>%
              dplyr::filter(site %in% nd)) %>%
  dplyr::select(-X) %>%
  group_by(site) %>%
  arrange(site, dist) %>%
  mutate(type = str_sub(site, 1, 1),
         num = str_sub(site, 2, 2))

df <- df %>%
  group_by(site) %>%
  transmute(range = gam >= (sill - 0.0005),
            dist = dist) %>%
  dplyr::filter(range == TRUE) %>%
  dplyr::filter(dist == min(dist)) %>%
  transmute(range = dist) %>%
  right_join(df)


semi_all <- ggplot(data = df, 
                   aes(x = dist, 
                       y = gam,
                       colour = num)) + 
  geom_hline(aes(yintercept = sill,
                 color = num),
             linetype = "dashed",
             show.legend = FALSE) + 
  geom_vline(aes(xintercept = range,
                 color = num),
             linetype = "dotted",
             show.legend = FALSE) +
  scale_colour_viridis(discrete = TRUE) + 
  scale_y_continuous(limits = c(0, 0.015),
                     breaks = seq(0, 0.015, 0.003)) +
  geom_point(shape = 1) + theme_bw() + 
  theme(legend.position = c(0.6, 0.72),
        aspect.ratio = 1,
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.margin =margin(t = 0.12, b = 0.12, l = 0.1, r = 0.1,
                              unit='cm'),
        legend.background = element_rect(
          fill = "gray90",
          linetype = "solid",
          colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  ylab("Semivariance") +
  xlab("Distance (m)") +
  facet_wrap(~type, ncol=3)

ggsave(plot = semi_all, 
       filename = "all_semivariograms_facet.tiff",
       device = "tiff",
       width = 4, height = 3, 
       units = "in")
