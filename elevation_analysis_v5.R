# 
# Author: Jake Diamond
# Purpose: Plot hummock heights from detrended data as a 
# function of median water table, noting which sites 
# have significant hummock differences
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load libraries
library(tidyverse)
library(broom)

# Load cleaned delineated hummock elevation data
elev <- read.csv("relative_elevations_all_v6.csv",
                 stringsAsFactors = FALSE)
elev$X <- NULL
# Load depth to confining layer data
conf <- read.csv("Soils/depth_to_confining.csv",
               stringsAsFactors = FALSE) %>%
  fill(plot)
# Load hummock hollow metadata
huho <- read.csv("hu.ho.csv", 
                 stringsAsFactors = FALSE)

# Data cleaning
huho$point <- paste(huho$plot, huho$position, sep = ".")
conf$point <- paste(conf$plot, conf$point, sep = ".")
huho[huho$site == "B1", "site"] <- "L1"
huho[huho$site == "B3", "site"] <- "L2"
huho[huho$site == "B6", "site"] <- "L3"
conf[conf$site == "B1", "site"] <- "L1"
conf[conf$site == "B3", "site"] <- "L2"
conf[conf$site == "B6", "site"] <- "L3"

# Join data
df <- left_join(huho, conf) %>%
  right_join(dplyr::filter(elev, 
                           !is.na(point)),
             by = c("site", "point", "plot"))

# Write data to R data file
# saveRDS(df, "elevations_all")

# Confining layer depth analysis ------------------------------------------
# Make sure that depth to confining layer is numeric
# and account for data >120 cm
df$depth <- ifelse(df$depth == ">120",
                   "150",
                   df$depth)
df$depth <- as.numeric(df$depth)
df$depth <- df$depth / 100

lm_model <- function(data) {
  lm(depth ~ value, data = data)
}

# First get data into long format
df_l <- df %>%
  dplyr::filter(!is.na(depth)) %>%
  dplyr::select(site, z, z_mod_lin, z_mod_quad, depth) %>%
  gather(key = "z_type", 
         value = "value", 
         -site, -depth)
# Linear analyses of confining depth vs elevation
lin_mods <- df_l %>%
  group_by(site, z_type) %>%
  nest() %>%
  mutate(model = map(data, lm_model),
         glance_lm = map(model, glance),  
         rsq = map_dbl(glance_lm, "r.squared"),
         pval = map_dbl(glance_lm, "p.value"),
         tidied = map(model, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

# Get x and y data for text on final plot
lin_mods <- lin_mods %>%
  dplyr::filter(z_type == "z") %>%
  mutate(x = c(0.36, 0.5, 0.3, 0.13, 0.0, 0.0, -0.2, 1.25, 0.35, 0.25),
         y = c(0.05, 0.3, 0.3, 0.6, 0.3, 0.18, 0.2, 1, 0.8, 0.75)) %>%
  right_join(lin_mods)

lin_mods$pval_text = ifelse(lin_mods$pval < 0.001, 
                                   "<0.001", 
                                   round(lin_mods$pval, 3))

# Write models to file
write.csv(lin_mods, "soil_versus_confining.csv")

# Plot confining depth vs soil elevation
p_conf <- ggplot(data = dplyr::filter(df_l,
                                      z_type == "z"),
                aes(x = depth,
                     y = value)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE) +
  ylab(expression("Relative soil elevation (m)")) +
  xlab("Depth to confining layer (m)") +
  theme_bw() +
  geom_abline(slope = 1,
              intercept = 0,
              color = "dark grey",
              linetype = "dashed") + 
  facet_wrap(~site, scales = "free") +
  geom_text(data = dplyr::filter(lin_mods,
                                 z_type == "z"),
            aes(x = x,
                y = y,
                label = paste0("p = ", pval_text)),
            show.legend = FALSE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 z_type == "z"),
            aes(x = x,
                y = y - 0.15,
                label = paste("list(R^2 ==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE)
p_conf

ggsave(plot = p_conf, 
       filename = "soil_elevation_vs_confining_layer.tiff",
       device = "tiff",
       dpi = 300)

# Site level hummock elevation analysis ---------------------------------------
# t-tests for hummock-hollow elevation difference 
# on validation points
ttests_site_val <- df %>%
  dplyr::filter(!is.na(point),
                !is.na(hu.ho)) %>%
  dplyr::select(site, point, hu.ho, z) %>%
  group_by(site) %>%
  nest() %>%
  mutate(data = map(data, spread, hu.ho, z)) %>%
  mutate(
    ttest = map(data, ~ t.test(.x$hummock, .x$hollow)),
    tidied = map(ttest, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE)


# Dataframe for p-values and estimates
df_text <- ttests_site_val %>%
  dplyr::select(site, p.value, estimate)
df_text$p.value <- ifelse(df_text$p.value < 0.001, 
                          "<0.001", 
                          round(df_text$p.value, 3))
df_text$y <- c(1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.1, 1.7)
df_text$x <- "hollow"

# Quick plot of hummock heights
p_huho <- ggplot(data = df, 
                 aes(x = hu.ho,
                     y = z_d_lin)) +
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
                label = paste("list(Delta*z ==", 
                              round(estimate, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE) +
  ylab("Elevation above mean water table (m)") +
  xlab("")
p_huho
ggsave(plot = p_huho,
       filename = "hummock_hollow_elevation_diffs.tiff",
       device = "tiff")

# Site level hummock height versus water table ----------------------------
# Use the z coord at the well as our datum
elev2 <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(site = site,
            z_well = z) %>%
  right_join(elev) %>%
  mutate(z = z - z_well,
         zmin = zmin - z_well,
         zmax = zmax - z_well)

ggplot(data = dplyr::filter(elev2, area > 0.1),
       aes(x = zmaxn))+
  geom_histogram()+
  facet_wrap(~site)

ggplot(data = dplyr::filter(elev2, point == "well"),
       aes(x = site,
           y = mean))+
  geom_col()

# Define lowland sites
lsites <- c("L1", "L2", "L3")

# Lowland sites are artifically high because the highest point
# on each hummock is actually a tree trunk...need to fix
elev2 <- elev2 %>%
  mutate(z = ifelse(site %in% lsites,
                    z20,
                    z),
         zmax = ifelse(site %in% lsites,
                       z80,
                       z))

# Find average hummock height of delineated hummocks
# Want to subset hummocks that are > 0.1 m2 to 
# remove obvious outliers
hu_hts_avg <- elev2 %>%
  dplyr::filter(is.na(point), area > 0.1) %>%
  dplyr::select(site, z, zmin, zmax, zmaxn, zmeann, z80n) %>%
  transmute(site = site,
            ht_abs_mean = z,
            ht_abs = zmax,
            ht_rel = zmaxn,
            ht_rel_mean = zmeann,
            ht_rel_80 = z80n
            ) %>%
  gather(key = "z_type",
         value = "value",
         -site) %>%
  group_by(site, z_type) %>%
  summarize(avg_ht = mean(value),
            sd_ht = sd(value)) %>%
  left_join(dplyr::select(elev, site, mean,
                          median, quant_80, variation) %>%
              distinct())
# Linear model
lm_model_hums <- function(data) {
  lm(avg_ht ~ value, data = data)
}

# Look at all correlations of hydrometrics to hummock site level averages
corr_hts <- hu_hts_avg %>%
  dplyr::select(-sd_ht) %>%
  # dplyr::filter(!(site %in% c("L1", "L2", "L3", "D2", "D4"))) %>%
  gather(key = "hydromet", value = "value", 
         -site, -avg_ht, -z_type) %>%
  group_by(hydromet, z_type) %>%
  nest() %>%
  mutate(model = map(data, lm_model_hums),
         glance_lm = map(model, glance),  
         rsq = map_dbl(glance_lm, "r.squared"),
         pval = map_dbl(glance_lm, "p.value"),
         tidied = map(model, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

# P value text for plot
corr_hts$pval_text = ifelse(corr_hts$pval < 0.001, 
                            "<0.001", 
                            round(corr_hts$pval, 3))

# Levels for facet headings
facets <- c(
  `ht_abs_mean` = "Site level hummock height",
  `ht_abs` = "Site level datum",
  `ht_mean` = "Local mean hummock height",
  `ht_rel` = "Local hollow datum"
)

# Plot average hummock height (raw z coord) versus average water 
# table (raw z coord)
avg_hu_p <- ggplot(data = hu_hts_avg,
                     # dplyr::filter(hu_hts_avg,
                     #                    !(site %in% c("L1", "L2", "L3", "D2", "D4"))),
                   aes(x = mean,
                       y = avg_ht,
                       color = z_type)) +
  geom_point() + 
  # geom_smooth(data = dplyr::filter(hu_hts_avg,
  #                                  z_type %in% 
  #                                    c("ht_mean")),
  #             method = "lm",
  #             se = FALSE,
  #             show.legend = FALSE) +
  geom_errorbar(aes(x = mean,
                    ymin = avg_ht - sd_ht,
                    ymax = avg_ht + sd_ht),
                show.legend = FALSE) + 
  geom_errorbarh(aes(y = avg_ht,
                     xmin = mean - sqrt(variation),
                     xmax = mean + sqrt(variation)),
                 show.legend = FALSE) + 
  theme_bw() +
  # scale_color_manual(name = "Height metric",
  #                    breaks = c("ht_abs", "ht_rel"),
  #                    labels = c("Elevation", "Relative elevation"),
  #                    values = c("darkblue", "lightblue")) +
  theme(legend.position = "none",
        legend.background = element_rect(fill = "lightgray",
                                         size = 0.5, 
                                         linetype = "solid", 
                                         colour ="darkgrey")) +
  # geom_text(data = dplyr::filter(corr_hts,
  #                                z_type == "ht_mean",
  #                                hydromet == "mean"),
  #           aes(x = 0.2,
  #               y = 0,
  #               label = paste0("y = ", round(value, 2), "x",
  #                              "+", round(`(Intercept)`, 2))),
  #           show.legend = FALSE) +
  # geom_text(data = dplyr::filter(corr_hts,
  #                                z_type == "ht_mean",
  #                                hydromet == "mean"),
  #           aes(x = 0.2,
  #               y = -0.11,
  #               label = paste0("p = ", pval_text)),
  #           show.legend = FALSE) +
  # geom_text(data = dplyr::filter(corr_hts,
  #                                z_type == "ht_mean",
  #                                hydromet == "mean"),
  #           aes(x = 0.2,
  #               y = -0.05,
  #               label = paste("list(R^2 ==",
  #                             round(rsq, digits=2), ")")),
  #           show.legend = FALSE,
  #           parse = TRUE) +
  xlab("Mean daily water table (m)") + 
  ylab("Hummock height (m)") +
  facet_wrap(~z_type, scales = "free")#, labeller = as_labeller(facets))
avg_hu_p

ggsave(plot = avg_hu_p,
       filename = "mean_hummock_ht_vs_water_table_v2.tiff",
       device = "tiff",
       dpi = 300)
# Individual hummocks vs water table analysis ----------------------------------------
# Need to determine which sites need which kind of detrending
# # Note: T1 is riddled with hummocks, not looking bimodal
# D2, L3 maybe needs to be detrended (quad)
# D1, L2 (maybe), T2, T3 (maybe) needs to be detrended (linear)
# D3, D4, L3 does not need to be detrended
# 
# 

# Use the z coord at the well as our datum
elev2 <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(site = site,
            z_well = z) %>%
  right_join(elev) %>%
  mutate(z = z - z_well,
         zmin = zmin - z_well,
         zmax = zmax - z_well)

# Get hydro data relative to detrended z-coords
elev_det <- elev2 %>%
  dplyr::filter(point == "well") %>%
  dplyr::select(mean, median, quant_80, 
                quant_20, z, z_mod_quad, site) %>%
  transmute(mean = mean - z_mod_quad,
            median = median - z_mod_quad,
            quant_80 = quant_80 - z_mod_quad,
            quant_20 = quant_20 - z_mod_quad,
            site = site) %>%
  right_join(dplyr::select(elev, -mean, -median, 
                           -quant_80, -quant_20))

# First, for each site, want to plot individual
# hummock height
# as a function of (base/hollow) distance from water table
# Do this assuming that the linear detrended surface is the hollow
df_hu <- elev_det %>%
  dplyr::filter(is.na(point),
                between(area, 0.1, 10)) %>%
  mutate(dist_hydmedian = (zmin - z_mod_quad) - median,
         type = str_sub(site, 1, 1),
         number = str_sub(site, 2, 2))
# Linear model for hummock heights as a function of distance from water table
lm_hu_hts <- function(data) {
  lm(zmeann ~ dist_hydmedian, data = data)
}

# Get linear relationships of hummock height vs water table distance
lms_hts <- df_hu %>%
  dplyr::select(site, type, number, dist_hydmedian, zmeann) %>%
  # gather(key = "hydromet", value = "value", 
  #        -site, -avg_ht, -z_type) %>%
  group_by(site, type, number) %>%
  nest() %>%
  mutate(model = map(data, lm_hu_hts),
         glance_lm = map(model, glance),  
         rsq = map_dbl(glance_lm, "r.squared"),
         pval = map_dbl(glance_lm, "p.value"),
         tidied = map(model, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

# Plot this
hums_p <- ggplot(data = df_hu,
                 aes(x = dist_hydmedian,
                     y = zmax,
                     color = number)) +
  geom_point() + 
  geom_smooth(method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  theme_bw() + 
  xlab("Distance from water table (m)") + 
  ylab("Hummock height (m)") + 
  # geom_text(data = lms_hts,
  #           aes(x = 0.4,
  #               y = 0.75,
  #               group = site,
  #               label = paste0("bar(m)==",
  #                              round(mean(dist_hydmedian), 2),
  #                              "%+-%",
  #                              round(sd(dist_hydmedian), 2))),
  #           parse = TRUE,
  #           show.legend = TRUE) +
  facet_wrap(~ site, scales = "free")

hums_p
