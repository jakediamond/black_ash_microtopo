# 
# Author: Jake Diamond
# Purpose: Plot all site semivariograms based on best detrends
# Date: February 5, 2019
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(viridis)
library(cowplot)
library(ggpubr)
library(tidyverse)
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
l <- c("L2")
q <- c("D1", "D4", "T1", "D2")

# Get all data
df <- read.csv("semi_1m_results_quad.csv",
               stringsAsFactors = FALSE) %>%
  dplyr::filter(site %in% q) %>%
  bind_rows(read.csv("semi_1m_results_linear.csv",
                     stringsAsFactors = FALSE) %>%
              dplyr::filter(site %in% l)) %>%
  bind_rows(read.csv("semi_1m_results_no_trend.csv",
                     stringsAsFactor = FALSE) %>%
              dplyr::filter(site %in% nd)) %>%
  dplyr::select(-X) %>%
  group_by(site) %>%
  arrange(site, dist) %>%
  mutate(type = str_sub(site, 1, 1),
         num = str_sub(site, 2, 2))

# Get ranges and sills for data (modeled)
sr <- read.csv("semi_fit_results_quad.csv",
                     stringsAsFactors = FALSE) %>%
  dplyr::filter(site %in% q) %>%
  bind_rows(read.csv("semi_fit_results_linear.csv",
                     stringsAsFactors = FALSE) %>%
              dplyr::filter(site %in% l)) %>%
  bind_rows(read.csv("semi_fit_results_no_trend.csv",
                     stringsAsFactor = FALSE) %>%
              dplyr::filter(site %in% nd)) %>%
  dplyr::select(-X) %>%
  arrange(site) %>%
  mutate(type = str_sub(site, 1, 1),
         num = str_sub(site, 2, 2))
# Calculate more realistic ranges because models arent that good
sr <- left_join(df, sr) %>%
  group_by(site) %>%
  transmute(range = gam >= (psill),
            dist = dist) %>%
  dplyr::filter(range == TRUE) %>%
  dplyr::filter(dist == min(dist)) %>%
  transmute(range2 = dist) %>%
  right_join(sr) %>%
  mutate(rangeuse = ifelse(site %in% c("L1", "L2", "L3") |
                             is.na(range2),
                           range,
                           range2))

# color blind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73")
# Plot
semi_all <- ggplot(data = df, 
                   aes(x = dist, 
                       y = gam,
                       colour = num)) + 
  geom_hline(data = sr,
             aes(yintercept = psill,
                 color = num,
                 linetype = "psill"),
             show.legend = TRUE) + 
  geom_segment(data = sr,
             aes(x = rangeuse,
                 xend = rangeuse,
                 y = 0,
                 yend = psill,
                 color = num,
                 linetype = "range"),
             show.legend = TRUE) +
  scale_colour_manual(values = cbbPalette) + 
  scale_y_continuous(limits = c(0, 0.015),
                     breaks = seq(0, 0.015, 0.003)) +
  scale_x_continuous(breaks = seq(0, 20, 5)) + 
  scale_linetype_manual(name = "Variogram parameter",
                        values = c("psill" = "dashed",
                                  "range" = "dotted"), 
                        labels = c("Partial sill",
                                   "Range")) +
  geom_point(size = 1.5) + 
  theme_bw() + 
  theme(legend.position = c(0.5, 0.73),
        aspect.ratio = 1,
        legend.key.size = unit(0.2, "cm"),
        legend.box = "horizontal",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.margin = margin(t = 0.12, b = 0.12, l = 0.1, r = 0.1,
                              unit='cm'),
        legend.background = element_rect(
          fill = "gray90",
          linetype = "solid",
          colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))
        ) +
  ylab("Semivariance") +
  xlab("Distance (m)") +
  facet_wrap(~type, ncol=3, scales = "free_x")
semi_all
ggsave(plot = semi_all, 
       filename = "all_semivariograms_facet_new_v3.tiff",
       device = "tiff",
       dpi = 300,
       width = 6, height = 4, 
       units = "in")
