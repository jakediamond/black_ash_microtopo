# 
# Author: Jake Diamond
# Purpose: To get elevations relative to water table for all points
# Date: September 3, 2018
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")

# Load Libraries
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(lubridate)
library(purrr)
library(stringr)

# Load raw elevation data
elev <- read.csv("Lidar/Point_validation/valpts_r.csv")

# Get all elevation data relative to each site's well
elev <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(z_well = z, site = site) %>%
  right_join(elev) %>%
  mutate(z_rel = z - z_well)

# Rename "b" sites to "l", make upper-case
elev$site <- toupper(elev$site)
elev[elev$site == "B1", "site"] <- "L1"
elev[elev$site == "B3", "site"] <- "L2"
elev[elev$site == "B6", "site"] <- "L3"

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
            site = site)

# Calculate relative elevation of points to average water table
elev_rel <- elev %>%
  left_join(hydro_well) %>%
  transmute(z_relh_mean = z_rel - mean,
            z_relh_median = z_rel - median,
            z_relh_80 = z_rel - quant_80,
            z_relh_20 = z_rel - quant_20,
            site = site,
            point_id = point_id,
            x = x,
            y = y,
            plot = plot,
            point = point)

# Save data
write.csv(elev_rel, "relative_elevations.csv")

