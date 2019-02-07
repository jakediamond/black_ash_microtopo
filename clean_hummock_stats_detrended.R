# 
# Author: Jake Diamond
# Purpose: Hummock statistics for cleaned and detrended z data
# Date: February 5, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(rgeos)
library(sp)


# Read in all files and merge
df <- read_rds("hummocks_clean_detrended")

# Hummock statistics function
humstats_fun <- function (data, res){
  return(data.frame(zmin = min(data$z),
                    zmax = max(data$z),
                    zmean = mean(data$z),
                    zmaxn = max(data$z - min(data$z)),
                    z20 = quantile(data$z, 0.2),
                    z80 = quantile(data$z, 0.8),
                    z20n = quantile(data$z - min(data$z), 0.2),
                    z80n = quantile(data$z - min(data$z), 0.8),
                    zmeann = mean(data$z - min(data$z))
                    )
         )
} 

# analysis of hummock stats
df_stat <- df %>%
  group_by(site, id, trend) %>%
  nest() %>%
  mutate(stats = map(data, humstats_fun, res = 0.01)) %>%
  unnest(stats, .drop = TRUE)

# Write data to disc
write.csv(df_stat, "hummock_stats_clean_detrended_z.csv")
write_rds(df_stat, "hummock_stats_clean_detrended_z")
