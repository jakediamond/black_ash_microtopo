# 
# Author: Jake Diamond
# Purpose: Combine all soil chem data (CN + extractions)
# Date: February 8, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)

# Load data
df <- read.csv("Soils/Soil Extraction Chemistry/meta_for_r.csv",
               stringsAsFactors = FALSE)
df_chem <- read.csv("Soils/Soil Extraction Chemistry/soil_chem_r.csv",
                    stringsAsFactors = FALSE)
df_cn <- read.csv("Soils/CN_r.csv")

# Load cleaned delineated hummock elevation data
elev <- read.csv("relative_elevations_all_v6.csv",
                 stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(point)) %>%
  dplyr::select(-X, -(id:ag_ht.cm),
                -(mean:ET.P_avg))

# Plot all data based on best detrend
trend <- data.frame(site = c("D1", "D2", "D3", "D4", 
                             "L1", "L2", "L3",
                             "T1", "T2", "T3"),
                    trend = c("q", "l", "nd", "q",
                              "nd", "l", "nd",
                              "q", "nd", "nd"))
elev <- left_join(elev, trend)

elev <- elev %>%
  mutate(zd = ifelse(trend == "q",
                     z - z_mod_quad,
                     ifelse(trend == "l",
                            z - z_mod_lin,
                            z)))
# Well adjustments based on detrending
well <- dplyr::filter(elev, point == "well") %>%
  dplyr::select(site, trend, z, zd) %>%
  mutate(well_adj = ifelse(trend == "nd", 0, zd))

# Get hydro data relative to detrended z-coords
hydro <- read.csv("average_wt_info_new_sites_v4.csv") %>%
  dplyr::select(-X, -(8:14)) %>%
  dplyr::filter(!(site %in% c("S1", "S2", "L4"))) %>%
  right_join(dplyr::select(well, site, trend, well_adj)) %>%
  distinct() %>%
  dplyr::select(-variation, -avghydro_rel, -avghydro, -ET.P_avg) %>%
  gather(key = "hydromet", value = "value",
         -site, -trend, -well_adj) %>%
  mutate(value_adj = value + well_adj) %>%
  dplyr::select(-value) %>%
  spread(key = hydromet,
         value = value_adj)

# Use the z coord at each site's well as our datum and 
# combine with hydro
elev <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(site = site,
            z_well = zd) %>%
  right_join(dplyr::filter(elev,
                           !is.na(point))) %>%
  mutate(zd = zd - z_well) %>%
  left_join(hydro)

# Get depth to confining layer
df <- read.csv("soils/depth_to_confining.csv",
               stringsAsFactors = FALSE) %>%
  fill(plot) %>%
  right_join(df)

# Make sure that depth to confining layer is numeric
# and account for data >120 cm
df$depth <- ifelse(df$depth == ">120",
                   NA,
                   df$depth)
df$depth <- as.numeric(df$depth)

# Change B sites to L sites
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"

# Change point columns to match elev format (Plot.Point)
df$point <- paste(df$plot, df$point, sep = ".")

# CN analysis -------------------------------------------------------------
# Subset CN data for what we need
df_cn <- df_cn %>%
  dplyr::select(Name,
                X..N,
                X..C,
                CNRatio) %>%
  rename(N = X..N,
         C = X..C,
         CN = CNRatio)

# Get sites in correct format for CN data
df_cn$site <- str_sub(df_cn$Name,
                      start = 1,
                      end = 2)

# Get replicates in similar format to soil extraction data
df_cn$Name <- as.character(df_cn$Name)
df_cn$reps <- str_trim(unlist(strsplit(df_cn$Name, "R")))

# Change B sites to L sites
df_cn[df_cn$site == "B1", "site"] <- "L1"
df_cn[df_cn$site == "B3", "site"] <- "L2"
df_cn[df_cn$site == "B6", "site"] <- "L3"

# Change point columns to match elev format (Plot.Point)
df_cn$point <- str_sub(df_cn$reps,
                       start = 4,
                       end = -1L)

# Find average and rpd of replicates
df_cn_rep <- df_cn %>%
  group_by(reps) %>%
  dplyr::filter(n() > 1) %>%
  group_by(reps) %>%
  summarize_at(vars(N:CN),
               funs(mean, diff)) %>%
  gather(key, value, -reps) %>%
  extract(key, c("solute", "calc"), "(^[^_]+)_(.*)$") %>%
  spread(calc, value) %>%
  mutate(rpd = abs(diff) * 100 / mean) %>%
  left_join(dplyr::select(df_cn,
                          site,
                          reps,
                          point) %>%
              distinct(site, point,
                       .keep_all = TRUE))

# Graph of rpd
# p_rep_cn <- ggplot(data = df_cn_rep, 
#                 aes(x = solute,
#                     y = rpd)) +
#   geom_bar(stat = "summary", 
#            fun.y = "mean",
#            position = "dodge") + 
#   stat_summary(fun.data = mean_se, 
#                geom = "errorbar",
#                width = 0.4,
#                position = "dodge") +
#   theme_bw() + 
#   facet_wrap(~site) +
#   xlab("Solute") + 
#   ylab("Relative percent difference (%)")

# p_rep_cn
# Save plot
# ggsave(plot = p_rep_cn, 
# "rpd_for_replicates_cn_site.tiff",
# device = "tiff",
# width = 8,
# height = 6,
# units = "in")

# Get all data in one place
df_cn_all <- df_cn %>%
  group_by(reps) %>%
  filter(n() == 1) %>%
  gather(solute, mean, -reps, -Name, -site, -point) %>%
  ungroup() %>%
  full_join(df_cn_rep) %>%
  select(-rpd, -diff, -Name, -reps)

# Add CN data to metadata
df_cn_spread <- df_cn_all %>%
  spread(solute, mean) %>%
  right_join(df)

# Put data in long format with metadata
df_cn_l <- df_cn_spread %>%
  gather(key = solute, 
         value = conc,
         C, N, CN)

# Soil Extraction analysis -------------------------------------------------
# Same thing for chemistry data
df_chem$sample_id <- as.character(df_chem$sample_id)
df_chem$reps <- unlist(strsplit(df_chem$sample_id, "R"))

# Make calcium numeric because it thinks its a character
df_chem$ca <- as.numeric(df_chem$ca)

# Find average and rpd of replicates
df_chem_rep <- df_chem %>%
  group_by(reps) %>%
  dplyr::filter(n() > 1) %>%
  group_by(reps) %>%
  summarize_at(vars(cl:po4),
               funs(mean, diff)) %>%
  gather(key, value, -reps) %>%
  extract(key, c("solute", "calc"), "(^[^_]+)_(.*)$") %>%
  spread(calc, value) %>%
  mutate(rpd = abs(diff) * 100 / mean) %>%
  mutate(reps = as.numeric(reps))

df_chem_rep2 <- df %>%
  rename(reps = sample_id) %>%
  right_join(df_chem_rep, by = "reps")

# Graph of rpd
# p_rep_chem <- ggplot(data = df_chem_rep2, 
#                 aes(x = solute,
#                     y = rpd)) +
#   geom_bar(stat = "summary", 
#            fun.y = "mean",
#            position = "dodge") + 
#   stat_summary(fun.data = mean_se, 
#                geom = "errorbar",
#                width = 0.4,
#                position = "dodge") +
#   theme_bw() + 
#   facet_wrap(~site) +
#   xlab("Solute") + 
#   ylab("Relative percent difference (%)")

# p_rep_chem
# Save plot
# ggsave(plot = p_rep_chem, 
# "rpd_for_replicates_chem_site.tiff",
# device = "tiff",
# width = 8,
# height = 6,
# units = "in")


# Get all data in one place
df_chem_all <- df_chem %>%
  mutate(reps = as.numeric(reps)) %>%
  group_by(reps) %>%
  filter(n() == 1) %>%
  gather(solute, mean, -reps, -sample_id) %>%
  full_join(df_chem_rep) %>%
  select(-rpd, -diff) %>%
  mutate(sample_id = reps,
         sample_id = as.integer(sample_id))

# Add CN data to metadata
df_chem_spread <- df_chem_all %>%
  spread(solute, mean) %>%
  right_join(df)

# Put data in long format with metadata
df_chem_l <- df_chem_spread %>%
  gather(key = solute, 
         value = conc,
         ca, cl, mg, no3, po4, so4)

# Overall soil analysis ---------------------------------------------------
# Combine CN and soil chem extraction data (long and wide)
df_data_l <- df_chem_l %>%
  ungroup() %>%
  dplyr::select(-reps) %>%
  bind_rows(ungroup(df_cn_l)) %>%
  left_join(elev, by = c("site", "plot", "point"))

df_data_w <- df_chem_spread %>%
  ungroup() %>%
  dplyr::select(-reps) %>%
  left_join(df_cn_spread) %>%
  left_join(elev, by = c("site", "plot", "point"))


# Write both to disc
write_rds(df_data_l, "soil_chem_long")
write_rds(df_data_w, "soil_chem_wide")
