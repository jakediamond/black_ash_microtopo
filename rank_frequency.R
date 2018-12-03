# 
# Author: Jake Diamond
# Purpose: To plot D2 hummock volumes
# Date: May 19, 2018
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/Spatial_data")

# Load Libraries
library(dplyr)
library(ggplot2)
library(scales)
library(broom)

# Read in data
df <- read.csv("val_ds.csv")

# Calculate rank frequency
df <- df %>%
  mutate(rank = cume_dist(desc(V_val_m3)))

# Fit exponential line
mod <- tidy(lm(rank ~ exp(V_val_m3), data = df))

# Plot
p_exp <- ggplot(data = df,
       aes(x = V_val_m3,
           y = rank)) +
  geom_point(size = 3,
             shape = 1) +
  stat_smooth(formula = y ~ exp(x)) +
  annotate("text",
           x = 0.001,
           y = 0.1,
           label = "p = 1.8E-18",
           size = 6) + 
  theme_bw() +
  scale_x_log10(
    labels = trans_format('log10', math_format(10 ^ .x)),
    limits = c(10 ^ -3.3, 10 ^ -0.3),
    breaks = c(10 ^ -3, 10 ^ -2, 10 ^ -1)) +
  scale_y_log10(
    labels = trans_format('log10', math_format(10 ^ .x)),
    limits = c(10 ^ -2, 10 ^ 0),
    breaks = c(10 ^ -2, 10 ^ -1, 10 ^ 0)) +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 18, colour = "black"),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.12, 0.17),
    plot.title = element_text(lineheight = .8, face = "bold"),
    legend.text = element_text(size = 16),
    legend.background = element_rect(
      fill = "gray90",
      size = 1,
      linetype = "solid",
      colour = "black")) +
  xlab(expression("Hummock volume, v (" * m ^ 3 * ")")) +
  ylab(expression("P (" * V>=v* ")"))

ggsave(p_exp, filename = "D2_hummock_exponential.tiff",
       device = "tiff",
       dpi = 600)  

