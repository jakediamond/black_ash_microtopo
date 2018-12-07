# 
# Author: Jake Diamond
# Purpose: To analyze soil chemistry data
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/Soils/Soil Extraction Chemistry")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data/Soils/Soil Extraction Chemistry")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)
library(viridis)

# Load data
df <- read.csv("meta_for_r.csv", stringsAsFactors = FALSE)
df2 <- read.csv("soil_chem_r.csv")
elev <- read.csv("C:/Users/diamo/Dropbox/Projects/EAB/Data/relative_elevations.csv")
elev$X <- NULL

# Change point column in df to match elev
df$point <- paste(df$plot, df$point, sep = ".")

# Change B sites to L sites
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"

# Combine elevation and meta data
df <- left_join(df, elev, by = c("site", "plot", "point"))

# Quick plot of hummock heights
ggplot(data = df,
       aes(x = hu.ho,
           y = z_relh_mean)) +
  geom_boxplot() + 
  facet_wrap(~site)


# Find replicates
df2$sample_id <- as.character(df2$sample_id)
df2$reps <- unlist(strsplit(df2$sample_id, "R"))

# Make calcium numeric because it thinks its a character
df2$ca <- as.numeric(df2$ca)

# Find average and rpd of replicates
df_rep <- df2 %>%
  group_by(reps) %>%
  filter(n() > 1) %>%
  group_by(reps) %>%
  summarize_at(vars(cl:po4),
               funs(mean, diff)) %>%
  gather(key, value, -reps) %>%
  extract(key, c("solute", "calc"), "(^[^_]+)_(.*)$") %>%
  spread(calc, value) %>%
  mutate(rpd = abs(diff) * 100 / mean) %>%
  mutate(reps = as.numeric(reps))

df_rep2 <- df %>%
  rename(reps = sample_id) %>%
  right_join(df_rep)

# Graph of rpd
p_rep <- ggplot(data = df_rep2, 
                aes(x = solute,
                    y = rpd,
                    fill = hu.ho)) +
  geom_bar(stat = "summary", 
           fun.y = "mean",
           position = "dodge") + 
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               width = 0.4,
               position = "dodge") +
  theme_bw() + 
  facet_wrap(~site) +
  xlab("Solute") + 
  ylab("Relative percent difference (%)")

p_rep
# Save plot
ggsave(plot = p_rep, 
       "rpd_for_replicates_huho_site.tiff",
       device = "tiff",
       width = 8,
       height = 6,
       units = "in")

# Get all data in one place
df3 <- df2 %>%
  mutate(reps = as.numeric(reps)) %>%
  group_by(reps) %>%
  filter(n() == 1) %>%
  gather(solute, mean, -reps, -sample_id) %>%
  full_join(df_rep) %>%
  select(-rpd, -diff) %>%
  mutate(sample_id = reps,
         sample_id = as.integer(sample_id))

df <- df3 %>%
  spread(solute, mean) %>%
  right_join(df)

df_l <- df %>%
  gather(key = solute, 
         value = conc,
         ca, cl, mg, no3, po4, so4)

# Strip panel names for plot
levels(df_l$solute) <- c("Ca^{`2+`}",
                         "Cl^{`-`}",
                         "Mg^{`2+`}",
                         "NO[3]^{`-`}",
                         "PO[4]^{`3-`}",
                         "SO[4]^{`2-`}")

# Plot
chem_p <- ggplot(data = dplyr::filter(df_l,
                                      # !(site %in% 
                                      #     c("L1", "L2", "L3")
                                      #   ),
                                      # site == "D1",
                                      !(solute == "cl" &
                                          conc > 30),
                                      log == 0
                                      ),
                 aes(x = z_rel.d - (median - (z_mod)),
                     y = conc,
                     colour = site)) +
  geom_point(aes(shape = hu.ho)) +
  geom_smooth() +
  scale_shape_manual(name = "",
                     values = c(1, 16)) +
  scale_colour_viridis(discrete = TRUE) +
  ylab(expression("Soil Extraction Concentration (mg "*L^{-1}*")")) +
  xlab("Relative elevation (m)") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~solute, scales = "free",
             labeller = label_parsed)
chem_p

chem_p2 <- ggplot(data = dplyr::filter(df_l, solute == "so4"),
                 aes(x = z_relh_median,
                     y = conc)) +
  geom_point() +
  ylab(expression("Soil Extraction Concentration (mg "*L^{-1}*")")) +
  xlab("Relative elevation (m)") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~site, scales = "free")
chem_p2

# t-tests
ttests_site <- df_l %>%
  select(site, solute, hu.ho, conc) %>%
  filter(hu.ho != "#N/A") %>%
  group_by(site, solute) %>%
  spread(key = hu.ho,
         value = conc,
         fill = NA) %>%
  group_by(site, solute) %>%
  do(tidy(t.test(.$hollow, .$hummock)))

ttests_globe <- df_l %>%
  select(site, solute, hu.ho, conc) %>%
  filter(hu.ho != "#N/A") %>%
  group_by(site, solute) %>%
  spread(key = hu.ho,
         value = conc,
         fill = NA) %>%
  group_by(solute) %>%
  do(tidy(t.test(.$hollow, .$hummock)))

ttests_globe2 <- df_l %>%
  filter(hu.ho != "#N/A",
         intermediate == 0,
         log == 0,
         edge == 0) %>%
  select(site, solute, hu.ho, conc) %>%
  group_by(site, solute) %>%
  spread(key = hu.ho,
         value = conc,
         fill = NA) %>%
  group_by(solute) %>%
  do(tidy(t.test(.$hollow, .$hummock)))



# Dataframe for p-values
df_p <- ttests_globe2 %>%
  select(solute, p.value) %>%
  mutate(hu.ho = ifelse(solute == "SO[4]^{`2-`}", "hummock", "hollow"))
df_p$y <- c(87.5, 2, 7, 0.8, 3.5, 10)

# Quick global plot
p_globe <- ggplot(data = filter(df_l, 
                                hu.ho != "#N/A",
                                intermediate == 0,
                                log == 0,
                                edge == 0), 
                  aes(x = hu.ho, 
                      y = conc,
                      fill = hu.ho)) + 
  geom_bar(stat = "summary", 
           fun.y = "mean") + 
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               width = 0.4) + 
  geom_text(data = df_p,
            aes(x = hu.ho,
                y = y,
                label = paste0("p = ", round(p.value, 3))),
            show.legend = FALSE) +
  ylab(expression("Soil Extraction Concentration (mg "*L^{-1}*")")) +
  xlab("") +
  scale_fill_manual(values = c("black", "grey")) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~solute, scales = "free",
             labeller = label_parsed)

p_globe

# Save Plot
ggsave(filename = "soil_chem_global_no_logs_no_int.tiff",
       plot = p_globe,
       device = "tiff",
       dpi = 600)

# By site
ggplot(data = filter(df_l, hu.ho != "#N/A"), 
       aes(x = hu.ho, 
           y = conc,
           fill = site)) + 
  geom_bar(position = "dodge",
           stat = "summary", 
           fun.y = "mean") + 
  # stat_summary(fun.data = mean_se, 
  #              geom = "errorbar",
  #              width = 0.4) + 
  facet_wrap(~solute, scales = "free")


