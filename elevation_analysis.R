# 
# Author: Jake Diamond
# Purpose: Plot hummock heights from detrended data as a function of median water table, noting which sites have significant hummock differences
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load libraries
library(tidyverse)
library(broom)

# Load data
elev <- read.csv("relative_elevations.csv",
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
  select(site, point, hu.ho, z) %>%
  group_by(site) %>%
  nest() %>%
  mutate(data = map(data, spread, hu.ho, z)) %>%
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
           y = z_d)) +
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
p
ggsave(plot = p,
       filename = "hummock_hollow_elevation_diffs.tiff",
       device = "tiff",
       width = 8,
       height = 6)

saveRDS(df, "elevations_wt")

# Relative hydro data
df_rel <- read.csv("average_wt_info_new_sites_v4.csv")

# Plot hummock heights vs elevation
# get data
pdat <- df_text %>%
  select(site, estimate, p.value) %>%
  left_join(df_rel) %>%
  mutate(site_type = str_sub(site, 1 , 1))

lm_model <- function(data) {
    lm(estimate ~ value, data = data)
}

tidy_mod_fun <- function(model) {
  glance <- glance(model) %>%
    select(adj.r.squared)
}

# Look at every correlation
corr <- pdat %>%
  select(-X, -p.value, -site_type) %>%
  gather(key = "hydromet", value = "value", -site, -estimate) %>%
  group_by(hydromet) %>%
  nest() %>%
  mutate(model = map(data, lm_model),
         glance_lm = map(model, glance),  
         rsq = map_dbl(glance_lm, "r.squared"),
         pval = map_dbl(glance_lm, "p.value"),
         tidied = map(model, tidy)
         ) %>% 
  unnest(tidied, .drop = TRUE) %>%
  select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

# Plot
p_hu.wt <- ggplot(data = pdat,
                  aes(x = mean,
                      y = estimate)) +
  theme_bw() + 
  ylab("Mean hummock height (m)") + 
  xlab("Mean water table (m)") +
  geom_point() + 
  geom_smooth(method = "lm") +
  geom_text(data = filter(corr, hydromet == "mean"),
            aes(x = -0.3,
                y = 0.18,
                label = paste0("p = ", round(pval, 3))),
            show.legend = FALSE) +
  geom_text(data = filter(corr, hydromet == "mean"),
            aes(x = -0.3,
                y =  0.2,
                label = paste("list(R^2 ==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE)
p_hu.wt

ggsave(plot = p_hu.wt, filename = "hummock_height_hydro.tiff",
       device = "tiff",
       dpi = 300)

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
