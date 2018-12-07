# 
# Author: Jake Diamond
# Purpose: To get elevations relative to water table for all points
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/HydroData")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)


# Load raw elevation data
elev <- read.csv("valpts_detrend_r.csv")
elev$X <- NULL

# Calculate detrended elevations
# Need to account for the fact that points are on 1.2m stakes
elev <- elev %>%
  dplyr::filter(point != "well") %>%
  mutate(z = z - 1.2,
         z_d = z - z_mod) %>%
  full_join(., dplyr::filter(elev, point == "well")) %>%
  mutate(z_d = ifelse(point == "well",
                      z - z_mod,
                      z_d))

# Get all elevation data relative to each site's well
elev <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(z_well = z, 
            z_well.d = z_d,
            site = site) %>%
  right_join(elev) %>%
  mutate(z_rel = z - z_well,
         z_rel.d = z_d - z_well.d)

# Load hydrology data
hydro <- read.csv("HydroData/average_wt_info_new_sites_v2.csv")

# Load well height data
wells <- read.csv("HydroData/well_ht_aboveground.csv")

# Change hydrodata to be relative to well top, not ground surface (meters)
hydro_well <- hydro %>%
  left_join(wells) %>%
  transmute(mean = mean - ag_ht.cm / 100,
            median = median - ag_ht.cm / 100,
            quant_80 = quant_80 - ag_ht.cm / 100,
            quant_20 = quant_20 - ag_ht.cm / 100,
            site = site) %>%
  na.omit()

# Calculate relative elevation of points to average water table
elev_rel <- elev %>%
  left_join(hydro_well)

# Save data
write.csv(elev_rel, "relative_elevations.csv")

