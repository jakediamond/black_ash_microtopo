# 
# Author: Jake Diamond
# Purpose: Plot hummock heights from detrended data as a 
# function of median water table, noting which sites 
# have significant hummock differences
# Date: September 3, 2018
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load libraries
library(broom)
library(cowplot)
library(ggpubr)
library(viridis)
library(tidyverse)

# Load cleaned delineated hummock elevation data
elev <- read.csv("relative_elevations_all_v6.csv",
                 stringsAsFactors = FALSE)
elev$X <- NULL
# Use the z coord at each site's well as our datum
elev <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(site = site,
            z_well = z) %>%
  right_join(elev) %>%
  mutate(z = z - z_well,
         zmin = zmin - z_well,
         zmax = zmax - z_well,
         z80 = z80 - z_well)

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

# Confining layer depth analysis ------------------------------------------
# Join data
df <- left_join(huho, conf) %>%
  right_join(dplyr::filter(elev, 
                           !is.na(point)),
             by = c("site", "point", "plot"))

# Make sure that depth to confining layer is numeric
# and account for data >120 cm
df$depth <- ifelse(df$depth == ">120",
                   NA,
                   df$depth)
df$depth <- as.numeric(df$depth)
df$depth <- df$depth / 100

# Account for L3 drop off in plot 3
df$z <- ifelse(df$site == "L3" &
                 df$plot == 3,
               df$z + 0.30,
               df$z)

lm_model <- function(data) {
  data$z = data$z - data$depth
  lm(depth ~ z, data = data)
}

# First get data into long format
df_l <- df %>%
  dplyr::filter(!is.na(depth)) %>%
  dplyr::select(site, z, depth, hu.ho)

# Site types and numbers
df_l$type <- str_sub(df_l$site, 1, 1)
df_l$num <- str_sub(df_l$site, 2, 2)

# Linear analyses of confining depth vs elevation
lin_mods <- df_l %>%
  group_by(site, type, num) %>%
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
# Separated by hummock hollow
lin_mods_huho <- df_l %>%
  group_by(site, type, num, hu.ho) %>%
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

# Pvalue text for plot
lin_mods$pval_text = ifelse(lin_mods$pval < 0.001, 
                                   "p<0.001", 
                                   paste0("p=", round(lin_mods$pval, 3)))
lin_mods_huho$pval_text = ifelse(lin_mods_huho$pval < 0.001, 
                            "p<0.001", 
                            paste0("p=", round(lin_mods_huho$pval, 3)))

# Write models to file
write.csv(lin_mods, "soil_thick_versus_confining_depth_noNAs_L3change.csv")
write.csv(lin_mods, "soil_thick_versus_confining_depth_noNAs_huho_L3change.csv")
# Get sites that have significant relationships
sig_sites <- lin_mods %>%
  dplyr::filter(pval < 0.01) %>%
  dplyr::select(site) %>%
  pull()

# Plot confining depth vs soil elevation
p_conf_d <- ggplot(data = dplyr::filter(df_l,
                                      type == "D"),
                aes(x = z - depth,
                     y = depth)) +
  geom_point(aes(shape = hu.ho)) +
  scale_shape_manual(breaks = c("hollow", "hummock"),
                     labels = c("Hollow", "Hummock"),
                     values = c(1, 16),
                     name = "Microsite") +
  ylab(expression("Soil thickness (m)")) +
  xlab("Confining layer elevation (m)") +
  theme_bw() +
  theme(legend.justification = "top",
        # legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  geom_abline(aes(slope = -1,
                  intercept = 0,
                  linetype = "Smooth surface (-1:1)"),
              color = "dark grey",
              show.legend = FALSE) + 
  geom_hline(aes(yintercept = 0.2,
                 linetype = "Subsurface reflection"),
             color = "dark grey") +
  scale_linetype_manual(name = "Surface soil model",
                        breaks = c("Smooth surface (-1:1)",
                                   "Subsurface reflection"
                                   ),
                        values = c("dashed",
                                   "dotted")) + 
  scale_y_continuous(limits = c(0, 1.2),
                     breaks = seq(0, 1.5, 0.4)) +
  scale_x_continuous(limits = c(-1, 0.1),
                     breaks = seq(-1.5, 0, 0.5)) +
  geom_text(data = dplyr::filter(lin_mods,
                                 # site %in% sig_sites,
                                 type == "D"),
            aes(x = -0.25,
                y = 0.95,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 # site %in% sig_sites,
                                 type == "D"),
            aes(x = -0.25,
                y = 1.1,
                label = paste("list(slope==",
                              round(z, digits=1), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 type == "D"),
            aes(x = -0.25,
                y = 0.75,
                label = pval_text),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, scales = "free_x", ncol = 4)
p_conf_da <- p_conf_d + theme(legend.position = "none")

p_conf_l <- ggplot(data = dplyr::filter(df_l,
                                        type == "L"),
                   aes(x = z - depth,
                       y = depth)) +
  geom_point(aes(shape = hu.ho)) +
  scale_shape_manual(breaks = c("hollow", "hummock"),
                     values = c(1, 16)) +
  ylab(expression("Organic soil thickness (m)")) +
  xlab("Confining layer elevation (m)") +
  theme_bw() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  geom_abline(slope = -1,
              intercept = 0,
              color = "dark grey",
              linetype = "dashed") + 
  geom_hline(yintercept = 0.2,
             linetype = "dotted",
             color = "dark grey") +
  scale_y_continuous(limits = c(0, 1.2),
                     breaks = seq(0, 1.5, 0.4)) +
  scale_x_continuous(limits = c(-1, 0.1),
                     breaks = seq(-1.5, 0, 0.5)) +
  geom_text(data = dplyr::filter(lin_mods,
                                 # site %in% sig_sites,
                                 type == "L"),
            aes(x = -0.25,
                y = 0.95,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 # site %in% sig_sites,
                                 type == "L"),
            aes(x = -0.25,
                y = 1.1,
                label = paste("list(slope==",
                              round(z, digits=1), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 type == "L"),
            aes(x = -0.25,
                y = 0.75,
                label = pval_text),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, scales = "free_x", ncol = 4)

p_conf_t <- ggplot(data = dplyr::filter(df_l,
                                        type == "T"),
                   aes(x = z - depth,
                       y = depth)) +
  geom_point(aes(shape = hu.ho)) +
  scale_shape_manual(breaks = c("hollow", "hummock"),
                     values = c(1, 16)) +
  ylab(expression("Soil thickness (m)")) +
  xlab("Mineral layer elevation (m)") +
  theme_bw() +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  geom_abline(slope = -1,
              intercept = 0,
              color = "dark grey",
              linetype = "dashed") + 
  geom_hline(yintercept = 0.2,
             linetype = "dotted",
             color = "dark grey") +
  scale_y_continuous(limits = c(0, 1.2),
                     breaks = seq(0, 1.5, 0.4)) +
  scale_x_continuous(limits = c(-1, 0.1),
                     breaks = seq(-1.5, 0, 0.5)) +
  geom_text(data = dplyr::filter(lin_mods,
                                 site %in% sig_sites,
                                 type == "T"),
            aes(x = -0.25,
                y = 0.95,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 site %in% sig_sites,
                                 type == "T"),
            aes(x = -0.25,
                y = 1.1,
                label = paste("list(slope==",
                              round(z, digits=1), ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lin_mods,
                                 type == "T"),
            aes(x = -0.25,
                y = 0.8,
                label = pval_text),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, scales = "free_x", ncol = 4)

legend <- cowplot::get_legend(p_conf_d)

p_conf <- ggdraw() +
  draw_plot(p_conf_da + rremove("x.text") + rremove("x.title") +
              rremove("y.title"), 
            x = 0.035, y = 0.685, width = 0.96, height = 0.32) +
  draw_plot(p_conf_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.39, width = 0.77, height = 0.32) +
  draw_plot(p_conf_t + rremove("y.title"),
            x = 0.035, y = 0, width = 0.7365, height = 0.41) +
  draw_grob(legend,
            0.84, 0.1, .3/3.3, 0.5)
# p_conf
ggsave(plot = p_conf,
       filename = "depth_vs_elevation_confining_spread_Lchange_new.tiff",
       device = "tiff",
       width = 6.5, height = 4,
       units = "in")

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
# Define lowland sites
lsites <- c("L1", "L2", "L3")

# Lowland sites are artifically high because the highest point
# on each hummock is actually a tree trunk...need to fix
elev2 <- elev %>%
  mutate(zmaxn = ifelse(site %in% lsites,
                    zmeann,
                    zmaxn),
         zmeann = ifelse(site %in% lsites,
                       z20n,
                       zmeann),
         z80n = ifelse(site %in% lsites,
                         zmeann,
                         z80n))

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
  dplyr::summarize(avg_ht = mean(value),
            sd_ht = sd(value)) %>%
  left_join(dplyr::select(elev, site, mean,
                          median, quant_80, variation) %>%
              distinct())


# Correction for D2 water table for missing data
# hu_hts_avg[hu_hts_avg$site == "D2", "mean"] <- -0.1

# Linear model
lm_model_hums <- function(data) {
  lm(avg_ht ~ value, data = data)
}

# Look at all correlations of hydrometrics to hummock site level averages
corr_hts <- hu_hts_avg %>%
  dplyr::select(-sd_ht) %>%
  # dplyr::filter(!(site %in% lsites)) %>%
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

write.csv(corr_hts, "all_site_level_hummock_water_table_correlations.csv")
write.csv(corr_hts, "all_site_level_hummock_water_table_correlations_noL.csv")
# P value text for plot
corr_hts$pval_text = ifelse(corr_hts$pval < 0.001, 
                            "p<0.001", 
                            paste0("p=", round(corr_hts$pval, 3)))

hu_hts_p <- dplyr::filter(hu_hts_avg,
                          z_type == "ht_rel_80") %>%
  ungroup()

# Plot relative hummock height versus average water table
avg_hu_p <- ggplot(data = hu_hts_p,
                   aes(x = mean,
                       y = avg_ht)) +
  geom_point(aes(color = site)) + 
  geom_smooth(method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  geom_errorbar(aes(x = mean,
                    ymin = avg_ht - sd_ht,
                    ymax = avg_ht + sd_ht,
                    color = site),
                show.legend = FALSE) + 
  geom_errorbarh(aes(y = avg_ht,
                     xmin = mean - sqrt(variation),
                     xmax = mean + sqrt(variation),
                     color = site),
                 show.legend = FALSE) + 
  scale_color_viridis(discrete = TRUE,
                      name = "") +
  theme_bw() +
  theme(legend.position = c(0.1, 0.6),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        # legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.background = element_rect(
          fill = "gray90",
          linetype = "solid",
          colour = "black")) +
  geom_text(data = dplyr::filter(corr_hts,
                                 z_type == "ht_rel_80",
                                 hydromet == "mean"),
            aes(x = -0.5,
                y = 0.3,
                label = paste0("y = ", round(value, 2), "x",
                               "+", round(`(Intercept)`, 2))),
            show.legend = FALSE) +
  geom_text(data = dplyr::filter(corr_hts,
                                 z_type == "ht_rel_80",
                                 hydromet == "mean"),
            aes(x = -0.5,
                y = 0.24,
                label = pval_text),
            show.legend = FALSE) +
  geom_text(data = dplyr::filter(corr_hts,
                                 z_type == "ht_rel_80",
                                 hydromet == "mean"),
            aes(x = -0.5,
                y = 0.27,
                label = paste("list(R^2 ==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE) +
  xlab("Mean daily water table (m)") + 
  ylab("Mean hummock height (m)")
avg_hu_p

ggsave(plot = avg_hu_p,
       filename = "mean_hummock_ht_vs_water_table_D2_z80.tiff",
       device = "tiff",
       dpi = 300)
# Individual hummocks vs water table analysis ----------------------------------------
# Get all data based on best detrend
df_i <- read_rds("hummock_stats_clean_detrended_z") %>%
  dplyr::filter(zmaxn > 0.10) %>%
  mutate(zmaxn = ifelse(site %in% c("L1", "L2", "L3"),
                        zmeann,
                        zmaxn),
         zmeann = ifelse(site %in% c("L1", "L2", "L3"),
                         z20n,
                         zmeann),
         type = str_sub(site, 1, 1),
         num = str_sub(site, 2, 2))

# Get hydro data relative to detrended z-coords
hydro <- read.csv("average_wt_info_new_sites_v4.csv") %>%
  dplyr::select(-X, -(8:14)) %>%
  dplyr::filter(!(site %in% c("S1", "S2", "L4"))) %>%
  left_join(dplyr::select(df_i, site, trend)) %>%
  distinct()

# Join data together
df_i <- left_join(df_i, hydro)
# First, for each site, want to plot individual
# hummock height (actually 80% of max hummock height)
# as a function of (base/hollow) distance from water table
df_hu <- df_i %>%
  mutate(dist_hydmean = mean - zmin)

# Average hh-wt by site
avgh <- df_hu %>%
  mutate(hhwt = z80n - dist_hydmean) %>%
  group_by(site) %>%
  dplyr::summarize(mean = mean(hhwt, na.rm = T),
                   sd = sd(hhwt, na.rm = T))

# Plot all sites together
# First get linear model of all sites
lm_all <- lm(z80n ~ dist_hydmean, data = df_hu)
summary(lm_all)
lm_site <- lm(z80n ~ dist_hydmean*site, data = df_hu)
summary(lm_site)
anova(lm_all, lm_site)
# Same thing without L sites
lm_all2 <- lm(z80n ~ dist_hydmean, data = filter(df_hu, type == "D",
                                                 site != "D4"))
summary(lm_all2)
lm_site2 <- lm(z80n ~ dist_hydmean*site, data = filter(df_hu, type == "D",
                                                       site != "D4"))
summary(lm_site2)
anova(lm_site2, lm_all2)
plot(lm_all2)
# Plot this
hums_p_all <- ggplot(data = filter(df_hu, type == "D"),
                   aes(x = dist_hydmean,
                       y = z80n)) +
  geom_point(shape = 1,
             aes(color = site)) + 
  geom_smooth(
              method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  theme_bw() + 
  scale_color_viridis(discrete = TRUE) +
  theme(legend.justification = "top",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  xlab("Local hummock mean water table (m)") + 
  ylab("Hummock height (m)") + 
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1))

hums_p_all

# Linear model for hummock heights as a function of distance from water table
lm_hu_hts <- function(data) {
  lm(z80n ~ dist_hydmean, data = data)
}

# Get linear relationships of hummock height vs water table distance
lms_hts <- df_hu %>%
  dplyr::select(site, type, num, dist_hydmean, z80n) %>%
  group_by(site, type, num) %>%
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

lms_hts$pval_text = ifelse(lms_hts$pval < 0.001, 
                            "p<0.001",
                           paste0("p=",
                                  round(lms_hts$pval, 3)))

# Get x and y data for text on final plot
lms_hts <- lms_hts %>%
  mutate(x = c(-0.17, -0.17, -0.17, -0.17, -0.375, 
               -0.375, -0.375, -0.275, -0.275, -0.275),
         y = 0.48)

# Write models to file
write.csv(lms_hts, "individual_hummock_ht_vs_water_table_reverse_meanL_nodt.csv")
# Get sites that have significant relationships
sig_sites <- lms_hts %>%
  dplyr::filter(pval < 0.01) %>%
  dplyr::select(site) %>%
  pull()

# Plot this
hums_p_d <- ggplot(data = dplyr::filter(df_hu,
                                       type == "D"),
                 aes(x = dist_hydmean,
                     y = z80n)) +
  geom_point(shape = 1) + 
  geom_smooth(data = dplyr::filter(df_hu,
                                   type == "D",
                                   site %in% sig_sites),
              method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  theme_bw() + 
  theme(legend.justification = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  xlab("Local hummock mean water table (m)") + 
  ylab("Hummock height (m)") + 
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(-0.25, 0.15),
                     breaks = seq(-0.2, 0.1, 0.1)) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "D"),
            aes(x = x,
                y = y - 0.065,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "D"),
            aes(x = x,
                y = y,
                label = paste("list(slope==",
                              round(dist_hydmean, digits=1), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 type == "D"),
            aes(x = x + 0.01,
                y = y - 0.14,
                label = pval_text),
            show.legend = FALSE,
            size = 2.5) +
  facet_wrap(~ site, scales = "free_x", ncol = 4)

hums_p_l <- ggplot(data = dplyr::filter(df_hu,
                                        type == "L"),
                   aes(x = dist_hydmean,
                       y = zmeann)) +
  geom_point(shape = 1) + 
  geom_smooth(data = dplyr::filter(df_hu,
                                   type == "L",
                                   site %in% sig_sites),
              method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  theme_bw() + 
  theme(legend.justification = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  xlab("Local hummock mean water table (m)") + 
  ylab("Hummock height (m)") + 
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
    scale_x_continuous(limits = c(-0.47, -0.2),
                     breaks = seq(-0.4, -0.2, 0.1)) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "L"),
            aes(x = x,
                y = y - 0.065,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "L"),
            aes(x = x,
                y = y,
                label = paste("list(slope==",
                              round(dist_hydmean, digits=1), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 type == "L"),
            aes(x = x,
                y = y - 0.14,
                label = pval_text),
            show.legend = FALSE,
            size = 2.5) +
  facet_wrap(~ site, scales = "free_x", ncol = 4)

hums_p_t <- ggplot(data = dplyr::filter(df_hu,
                                        type == "T"),
                   aes(x = dist_hydmean,
                       y = z80n)) +
  geom_point(shape = 1) + 
  geom_smooth(data = dplyr::filter(df_hu,
                                   type == "T",
                                   site %in% sig_sites),
              method = "lm",
              se = FALSE,
              show.legend = FALSE) +
  theme_bw() + 
  theme(legend.justification = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  xlab("Local hummock mean water table (m)") + 
  ylab("Hummock height (m)") + 
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(-0.35, 0.07),
                     breaks = seq(-0.3, 0.0, 0.1)) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "T"),
            aes(x = x,
                y = y - 0.065,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 site %in% sig_sites,
                                 type == "T"),
            aes(x = x,
                y = y,
                label = paste("list(slope==",
                              round(dist_hydmean, digits=1), ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  geom_text(data = dplyr::filter(lms_hts,
                                 type == "T"),
            aes(x = x,
                y = y - 0.14,
                label = pval_text),
            show.legend = FALSE,
            size = 2.5) +
  facet_wrap(~ site, scales = "free_x", ncol = 4)

hums_p <- ggdraw() +
  draw_plot(hums_p_d + rremove("x.title") +
              rremove("y.title"), 
            x = 0.0347, y = 0.67, width = 0.96, height = 0.33) +
  draw_plot(hums_p_l + rremove("x.title"),
            x = 0, y = 0.36, width = 0.77, height = 0.33) +
  draw_plot(hums_p_t + rremove("y.title"),
            x = 0.0343, y = 0, width = 0.7365, height = 0.38)

ggsave(plot = hums_p,
       filename = "individ_hummock_vs_water_table_reverse_L_nodt.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")
