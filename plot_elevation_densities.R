# 
# Author: Jake Diamond
# Purpose: Plot all site density elevations based on best detrends
# Date: February 4, 2019
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
# No detrend sites
nd <- c("L3", "L1", "D3", "T2", "T3")
l <- c("L2", "D2")
q <- c("D1", "D4", "T1")

data_path <- "RDS_files/"
files <- dir(data_path, pattern = "*")
files <- files[!(files %in% "Samples")]

# Read in all files and merge
df <- tibble(filename = files,
             site = str_sub(filename, 1, 2),
             trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% nd,
                trend == "end") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "no")

df <- tibble(filename = files,
                site = str_sub(filename, 1, 2),
                trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% l,
                trend == "ear") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "linear") %>%
  bind_rows(df)

df <- tibble(filename = files,
             site = str_sub(filename, 1, 2),
             trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% q,
                trend == "uad") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "quad") %>%
  bind_rows(df)

# Write summary of this hummock data to disc for later
# df %>%
#   dplyr::filter(!is.na(id)) %>%
#   write_rds("hummocks_clean_detrended")

# Summarize data
df_sum <- read.csv("Lidar/hummock_stats_ext6.csv") %>%
  dplyr::select(id, site_area) %>%
  mutate(site = str_sub(id, 1, 2)) %>%
  dplyr::select(site, site_area) %>%
  distinct() %>%
  right_join(read.csv("hummock_stats_clean.csv",
                   stringsAsFactors = FALSE)) %>%
  dplyr::select(-X) %>%
  dplyr::filter(area > 0.1,
                zmaxn > 0.15) %>%
  mutate(zmaxn = ifelse(site %in% c("B1", "B3", "B6"),
                        z80n,
                        zmaxn)) %>%
  group_by(site) %>%
  summarize(no = n(),
            aratio = sum(area) / mean(site_area),
            zavg = mean(zmaxn),
            zsd = sd(zmaxn)) %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site))))

df$type <- str_sub(df$site, 1, 1)
df$num <- str_sub(df$site, 2, 2)
df_sum$type <- str_sub(df_sum$site, 1, 1)
df_sum$num <- str_sub(df_sum$site, 2, 2)

# Subsample data for faster processing
dfSam <- df %>%
  group_by(site) %>%
  sample_n(size = 10000)
rm(df)
# Plot data binary
p_hum_best_d <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "D"),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 8500),
                     breaks = seq(0, 8000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.43,
            y = 7900,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.53,
            y = 6850,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.48,
            y = 5700,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.justification = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_da <- p_hum_best_d + theme(legend.position = "none")

p_hum_best_l <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "L"),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 10000),
                     breaks = seq(0, 10000, 2500)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.43,
            y = 9400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.53,
            y = 8200,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.48,
            y = 6900,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
theme(legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      axis.text.y = element_blank()
) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_t <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "T"),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = seq(0, 8000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.43,
            y = 7400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.53,
            y = 6500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.48,
            y = 5450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

legend <- cowplot::get_legend(p_hum_best_d)

p_hum_best2 <- ggdraw() +
  draw_plot(p_hum_best_da + rremove("x.text") + rremove("x.title"), 
            x = 0, y = 0.7, width = 1, height = 0.3) +
  draw_plot(p_hum_best_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.4, width = 0.76, height = 0.3) +
  draw_plot(p_hum_best_t,
            x = 0, y = 0, width = 0.76, height = 0.4) + 
  draw_grob(legend,
            0.82, 0, .3/3.3, 0.5)

ggsave(plot = p_hum_best2,
       filename = "all_sites_densities_binary.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")

# Plot data with hollows classified
p_hum_best_d2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "D",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 7000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.justification = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_da2 <- p_hum_best_d2 + theme(legend.position = "none")

p_hum_best_l2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "L",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 6000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_t2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "T",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 7000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

legend2 <- cowplot::get_legend(p_hum_best_d)

p_hum_best3 <- ggdraw() +
  draw_plot(p_hum_best_da2 + rremove("x.text") + rremove("x.title"), 
            x = 0, y = 0.7, width = 1, height = 0.3) +
  draw_plot(p_hum_best_l2 + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.4, width = 0.76, height = 0.3) +
  draw_plot(p_hum_best_t2,
            x = 0, y = 0, width = 0.76, height = 0.4) + 
  draw_grob(legend2,
            0.82, 0, .3/3.3, 0.5)

ggsave(plot = p_hum_best3,
       filename = "all_sites_densities_hollow_classified.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")
