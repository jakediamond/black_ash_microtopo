# 
# Author: Jake Diamond
# Purpose: Plot hummock heights from detrended data as a 
# function of median water table, noting which sites 
# have significant hummock differences
# Date: February 3
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load libraries
library(tidyverse)
library(broom)

# Need to determine which sites need which kind of detrending
# # Note: T1 is riddled with hummocks, not looking bimodal
# D2, L3 maybe needs to be detrended (quad)
# D1, L2 (maybe), T2, T3 (maybe) needs to be detrended (linear)
# D3, D4, L3 does not need to be detrended
# 

no <- read.csv("mcl_results_no_trend.csv",
               stringsAsFactors = FALSE)
lin <- read.csv("mcl_results_linear.csv",
                stringsAsFactors = FALSE)
qua <- read.csv("mcl_results_quad.csv",
                stringsAsFactors = FALSE)

df <- no %>%
  dplyr::filter(site %in% 
                  c("D3", "D4", "L3")) %>%
  bind_rows(lin %>% 
              dplyr::filter(site %in% 
                              c("D1", "T2", "T3"))) %>%
  bind_rows(qua %>% 
              dplyr::filter(site %in% 
                              c("D2", "L1", "L2", "T1"))) %>%
  dplyr::select(G, mean, proportion, site, model)

df_dif <- df %>%
  group_by(site) %>%
  transmute(dif = max(mean, na.rm = TRUE) - 
              min(mean, na.rm = TRUE)) %>%
  distinct()


# Load hydrology data
hydro <- read.csv("average_wt_info_new_sites_v4.csv") %>%
  dplyr::select(-X, -(8:14))

df_dif <- left_join(df_dif, hydro)

ggplot(data = df_dif,
       aes(x = median,
           y = dif)) + geom_point()
