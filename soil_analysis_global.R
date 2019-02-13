# 
# Author: Jake Diamond
# Purpose: To analyze soil chemistry data globally
# Date: February 8, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/Soils/Soil Extraction Chemistry")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")
# Load Libraries
library(tidyverse)
library(broom)
library(ggplot2)
library(lubridate)
library(viridis)

# Load soil chemistry and data
df <- read_rds("soil_chem_long")

# Some outlier functions
# Z score
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}
# Median absolute deviation
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}
# Tukey's fences
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

# Apply functions
df <- df %>%
  group_by(site, solute) %>%
  mutate(insz = isnt_out_z(conc),
         insmad = isnt_out_mad(conc),
         instuk = isnt_out_tukey(conc),
         outlier = ifelse((insz + insmad +instuk) < 2,
                          1,
                          0))

# Calculate relative concentrations (relative to site average)
df <- df %>%
  filter(outlier == 0) %>%
  group_by(site, solute) %>%
  mutate(conc_avg = mean(conc, na.rm = TRUE),
         conc_rel = conc / conc_avg)

# Get a count of each measurement by site
df_count <- df %>%
  select(site, solute, conc, hu.ho) %>%
  group_by(site, hu.ho, solute) %>%
  summarize(cont = n())
write.csv(df_count, "site_level_soil_counts.csv")

# Get summary statistics for each site
df_sum <- df %>%
  select(site, solute, conc, sample_id) %>%
  group_by(site) %>%
  spread(solute, conc) %>%
  group_by(site) %>%
  summarise_all(funs(mean,sd), na.rm = T)
write.csv(df_sum, "site_level_soil_summary.csv")

# Get summary statistics for each site, split by hummock and hollow
df_sum2 <- df %>%
  select(site, solute, conc, sample_id, hu.ho) %>%
  group_by(site, hu.ho) %>%
  spread(solute, conc) %>%
  group_by(site, hu.ho) %>%
  summarise_all(funs(mean,sd), na.rm = T)
write.csv(df_sum2, "site_level_soil_summary_huho.csv")

# Define L sites
lsites <- c("L1", "L2", "L3")

# t-tests
# t-test based broken out by site
ttests_site <- df %>%
  filter(hu.ho != "#N/A") %>%
  select(site, solute, hu.ho, conc, sample_id) %>%
  group_by(site, solute) %>%
  spread(key = hu.ho,
         value = conc,
         fill = NA) %>%
  group_by(site, solute) %>%
  do(tidy(t.test(.$hollow, .$hummock)))
# Global test 
ttests_globe <- df %>%
  dplyr::filter(hu.ho != "#N/A") %>%
  dplyr::select(site, solute, hu.ho, conc_rel, sample_id) %>%
  group_by(site, solute) %>%
  spread(key = hu.ho,
         value = conc_rel,
         fill = NA) %>%
  group_by(solute) %>%
  do(tidy(t.test(.$hollow, .$hummock)))

# Dataframe for p-values
df_p <- ttests_globe %>%
  select(solute, p.value, estimate1, estimate2) %>%
  mutate(ratio = estimate1 / estimate2,
         hu.ho = "hollow",
         y = 0.4)

# Strip panel names for plot
df$solute <- factor(df$solute,
                           levels = c("ca",
                                      "cl",
                                      "mg",
                                      "no3",
                                      "po4",
                                      "so4",
                                      "C",
                                      "N",
                                      "CN"))
levels(df$solute) <- c("Ca^{`2+`}",
                              "Cl^{`-`}",
                              "Mg^{`2+`}",
                              "NO[3]^{`-`}",
                              "PO[4]^{`3-`}",
                              "SO[4]^{`2-`}",
                              "'%'*C",
                              "'%'*N",
                              "C:N")
df_p$solute <- factor(df_p$solute,
                    levels = c("ca",
                               "cl",
                               "mg",
                               "no3",
                               "po4",
                               "so4",
                               "C",
                               "N",
                               "CN"))
levels(df_p$solute) <- c("Ca^{`2+`}",
                       "Cl^{`-`}",
                       "Mg^{`2+`}",
                       "NO[3]^{`-`}",
                       "PO[4]^{`3-`}",
                       "SO[4]^{`2-`}",
                       "'%'*C",
                       "'%'*N",
                       "C:N")
# p value text
df_p$pvaltext <- ifelse(df_p$p.value < 0.001,
                        "p<0.001",
                        paste0("p=", round(df_p$p.value, 3)))

# Quick global plot
p_globe_rel <- ggplot(data = filter(df, 
                                hu.ho != "#N/A"), 
                  aes(x = hu.ho, 
                      y = conc_rel,
                      fill = hu.ho),
                  alpha = 0.8) + 
  geom_bar(stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               width = 0.4) + 
  geom_text(data = df_p,
            aes(x = hu.ho,
                y = y - 0.2,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = df_p,
            aes(x = hu.ho,
                y = y + 0.1,
                label = paste0("[Ho]:[Hu]=", 
                               round(ratio, 2))),
            show.legend = FALSE,
            size = 2) +
  ylab("Relative concentration (-)") +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  scale_fill_manual(breaks = c("hollow", "hummock"),
                    values = c("white", "dark grey")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))) +
  facet_wrap(~solute,
             labeller = label_parsed)

# Save Plot
ggsave(filename = "soil_chem_global_no_logs_no_int_rel.tiff",
       plot = p_globe_rel,
       device = "tiff",
       dpi = 300)

# By site
ggplot(data = filter(df, hu.ho != "#N/A"), 
       aes(x = hu.ho, 
           y = conc_rel,
           fill = site)) + 
  geom_bar(position = "dodge",
           stat = "summary", 
           fun.y = "mean") + 
  # stat_summary(fun.data = mean_se, 
  #              geom = "errorbar",
  #              width = 0.4) + 
  facet_wrap(~solute, scales = "free")


