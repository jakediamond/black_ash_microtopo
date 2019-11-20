# 
# Author: Jake Diamond
# Purpose: To analyze veg continuously with elev
# Date: December 6, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(vegan)

# Load veg data
df <- readRDS("veg_data.RDS")
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"

# Load elevation data
elev3 <- readRDS("elevations_wt")

# Diversity and richness analysis -----------------------------------------
# Get veg data into format for vegan diversity function...no moss
df_div <- df %>%
  filter(moss != 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Get moss data into format for vegan diversity function
df_moss <- df %>%
  filter(moss == 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Function to compute diversity
div_fun <- function(data) {
  badcols <- c("site", "hu.ho", "point")
  mat <- data %>%
    ungroup() %>%
    dplyr::select(-one_of(badcols))
  shann <- diversity(mat)
  simp <- diversity(mat, "simpson")
  rich <- sum(mat > 0)
  x <- data.frame(shannon = shann, 
                  simpson = simp, 
                  richness = rich)
  return(x)
}

# Calculate diversity by site and hummock/hollow
div_site <- df_div %>%
  group_by(site, hu.ho) %>%
  do(div_fun(.))

div_point <- df_div %>%
  group_by(site, point) %>%
  do(div_fun(.)) %>%
  left_join(df %>% 
              dplyr::select(site, point, depth) %>%
              unique()) %>%
  left_join(elev3, by = c("site", "point")) %>%
  ungroup()

div_point$depth2 <- as.numeric(as.character(div_point$depth))
div_point$edge <- ifelse(is.na(div_point$edge), 0, 1)
div_point$log <- ifelse(is.na(div_point$log), 0, 1)
div_point$notree <- ifelse(is.na(div_point$notree), 0, 1)
div_point$intermediate <- ifelse(is.na(div_point$intermediate), 
                                 0, 
                                 1)

# Model
lm_fun <- function(df) {
  lm(richness ~ z_relh_mean, data = df)
}

# Get linear model info
lms <- div_point %>%
  group_by(site) %>%
  nest() %>%
  mutate(model = map(data, lm_fun)) %>%
  mutate(
    adj.r.squared = map_dbl(model, 
                            ~ signif(summary(.)$adj.r.squared, 
                                     3)),
    intercept = map_dbl(model, ~ signif(.$coef[[1]], 3)),
    slope = map_dbl(model, ~ signif(.$coef[[2]], 3)),
    pvalue = map_dbl(model, ~ signif(summary(.)$coef[2, 4], 3)) 
    ) %>%
  dplyr::select(-data, -model)

colnames(lms) <- c("site", "rsq", "b", "m", "pval")
  
lms$ptext <- ifelse(lms$pval < 0.001, 
                          "<0.001", 
                          paste0("=", round(lms$pval, 3)))

# Quick plot
nomoss_elev <- ggplot(div_point, aes(x = z_relh_mean, 
                                     y = richness)) + 
  geom_point(shape = 21, 
             aes(fill = hu.ho,
                 alpha = factor(log))) +
  geom_smooth(method = "lm") + 
  geom_text(data = lms,
    aes(x = 1.4, y = 10, 
        label = paste0("list(y==~", m, "*x", b, ", ",
                       "~italic(R)^2==", rsq,
                       ")")),
    parse = TRUE) +
  geom_text(data = lms,
            aes(x = 1.4, y = 8,
                label = paste0("p ", ptext))) +
  scale_alpha_manual(breaks = c(0, 1),
                    values = c(1, 0.4)) +
  scale_fill_manual(breaks = c("ho", "hu"),
                      values = c("white", "black")) +
  theme_bw() + 
  xlab("Relative elevation (m)") +
  ylab("Richness") +
  theme(legend.position = c(0.75, 0.2)) + 
  facet_wrap(~site, scales = "free_x")

nomoss_elev

