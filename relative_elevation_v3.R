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

# Load raw elevation data for measured points
# elev <- read.csv("valpts_detrend_r.csv")
elev <- read.csv("clean_hummocks_and_valpts_detrend.csv") %>%
  dplyr::select(-X, -x, -xmin, -xmax, -y, -ymin, -ymax)

# Load well height data
wells <- read.csv("HydroData/well_ht_aboveground.csv")

# Load hydrology data
hydro <- read.csv("average_wt_info_new_sites_v4.csv") %>%
  dplyr::select(-X, -(8:14))

# Need to first account for the fact that points 
# are on 1.2m stakes, and also subtract well heights 
# from wells to get ground elevation at those 
elev <- elev %>%
  dplyr::filter(!is.na(point)) %>%
  left_join(wells) %>%
  mutate(z = ifelse(point == "well",
                    z - ag_ht.cm / 100,
                    z - 1.2)) %>%
  full_join(., dplyr::filter(elev, is.na(point)))

# Join data
elev_rel <- hydro %>%
  right_join(elev)

# Save data
# write.csv(elev_rel, "relative_elevations.csv")
# write.csv(elev_rel, "relative_elevations_quad.csv")
write.csv(elev_rel, "relative_elevations_all_v6.csv")
