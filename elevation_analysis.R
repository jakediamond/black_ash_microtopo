# 
# Author: Jake Diamond
# Purpose: Initial analysis of elevation data
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(lubridate)
library(purrr)
library(stringr)

# Load data
# elev <- read.csv("E:/Dropbox/Dropbox/Projects/EAB/Data/relative_elevations.csv")
elev <- read.csv("C:/Users/diamo/Dropbox/Projects/EAB/Data/relative_elevations.csv",
                 stringsAsFactors = FALSE)
elev$X <- NULL
huho <- read.csv("hu.ho.csv", stringsAsFactors = FALSE)

# Data cleaning
huho$point <- paste(huho$plot, huho$position, sep = ".")
huho[huho$site == "B1", "site"] <- "L1"
huho[huho$site == "B3", "site"] <- "L2"
huho[huho$site == "B6", "site"] <- "L3"

# Join data
df <- left_join(huho, elev, by = c("site", "plot" , "point"))

# t-tests for hummock-hollow elevation difference
ttests_site <- df %>%
  select(site, point, hu.ho, z_relh_mean) %>%
  group_by(site) %>%
  nest() %>%
  mutate(data = map(data, spread, hu.ho, z_relh_mean)) %>%
  mutate(
    ttest = map(data, ~ t.test(.x$hummock, .x$hollow)),
    tidied = map(ttest, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE)

# Dataframe for p-values and estimates
df_text <- ttests_site %>%
  select(site, p.value, estimate)
df_text$p.value <- ifelse(df_text$p.value < 0.001, 
                          "<0.001", 
                          round(df_text$p.value, 3))
df_text$y <- c(1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.1, 1.7)
df_text$x <- "hollow"

# Quick plot of hummock heights
p <- ggplot(data = df, 
       aes(x = hu.ho,
           y = z_relh_mean)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site) +
  geom_text(data = df_text,
            aes(x = x,
                y = y,
                label = paste0("p = ", p.value)),
            show.legend = FALSE) +
  geom_text(data = df_text,
            aes(x = x,
                y = y - 0.13,
                label = paste("list(Delta*z ==", round(estimate, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE) +
  ylab("Elevation above mean water table (m)") +
  xlab("")

ggsave(plot = p,
       filename = "hummock_hollow_elevation_diffs.tiff",
       device = "tiff",
       width = 8,
       height = 6)

saveRDS(df, "elevations_wt")

# Quick plot elev
depth <- ggplot(data = df, 
                aes(x = z_relh_mean, 
                    y = depth2 / 100)) + 
  geom_point() +
  theme_bw() + 
  xlab("Relative elevation (m)") +
  ylab("depth of peat (cm)") +
  theme(legend.position = "none") + 
  facet_wrap(~site) + 
  geom_abline(slope = 1, intercept = -1)

depth


div_point$depth3 <- ifelse(is.na(div_point$depth2), 130, div_point$depth2)
ggplot(div_point, aes(x = z_relh_mean, 
                      fill = plot, 
                      group = plot)) + 
  geom_density(alpha = 0.6) +
  facet_wrap(~site, scales = "free")
