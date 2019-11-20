# 
# Author: Jake Diamond
# Purpose: To analyze richness/diversity at the site level
# Date: February 27, 2019
# 

# Set Working Directory
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# Load Libraries
library(broom)
library(vegan)
library(Hmisc)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(tidyverse)
# Set plot theme
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))))
# Load veg data
df <- read_rds("veg_data.rds") %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site))))

# Get data into format for vegan diversity function, and calc richness
df_div_all <- df %>%
  dplyr::filter(moss == 0) %>%
  # mutate(tot = ifelse(is.na(number),
  #                     percent,
  #                     number)) %>%
  group_by(site, plot, species) %>%
  dplyr::summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Function to compute diversity
div_fun <- function(data) {
  badcols <- c("site")
  mat <- data %>%
    ungroup() %>%
    dplyr::select(-one_of(badcols))
  shann <- diversity(mat)
  simp <- diversity(mat, "simpson")
  rich <- sum(mat > 0)
  x <- data.frame(shannon = shann, simpson = simp, richness = rich)
  return(x)
}

# Caclulate diversity by site
div <- df_div_all %>%
  group_by(site, plot) %>%
  do(div_fun(.))

# Load forestry-hydro data
hyd <- read_csv("Forestry/forest_hydro_summary.csv")

# combine and plot
div <- left_join(div, hyd) %>%
  dplyr::select(median, site, plot, richness, shannon) %>%
  rename(Richness = richness, `Shannon diversity` = shannon) %>%
  gather(div, val, -site, -median, -plot)

# Linear model
lm_model <- function(data) {
  lm(depth ~ z, data = data)
}

# Linear analyses of confining depth vs elevation
lin_mods <- div %>%
  group_by(div) %>%
  nest() %>%
  mutate(model = map(data, ~lm(.$val ~ .$median)),
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
# Site types and numbers
div$type <- str_sub(div$site, 1, 1)
div$num <- str_sub(div$site, 2, 2)

# x and y data for lin mods text in graph
lin_mods$x <- 0.1
lin_mods$y1 <- c(27, 2.55)
lin_mods$y2 <- c(26, 2.5)
lin_mods$y3 <- c(24.5, 2.43)

# Plot
veghyd <- ggplot(data = div,
       aes(x = median,
           y = val)) + 
  stat_summary(fun.data = mean_cl_boot,
               aes(color = type)) + 
  geom_text_repel(aes(label = site),
                  nudge_x = 0.02,
                  stat = "summary",
                  fun.y=mean) +
  stat_smooth(method = "lm") + 
  scale_color_viridis(discrete = TRUE,
                      name = "Site type") +
  facet_wrap(~div, scales = "free_y") + 
  theme(legend.position = c(0.08, 0.20),
        legend.background = element_rect(fill = "grey90",
                                         color = "black"),
        legend.margin = margin(0.08,0.1,0.08,0.1, "cm")) +
  geom_text(data = lin_mods,
            aes(x = x,
                y = y2,
                label = paste("list(R^2==",
                              round(rsq, digits=2), ")")),
            show.legend = FALSE,
            size = 3,
            parse = TRUE) +
  geom_text(data = lin_mods,
            aes(x = x,
                y = y1,
                label = paste("list(y==",
                              round(`.$median`, digits=1),
                              "*x*+",
                              round(`(Intercept)`, digits = 1),
                              ")")),
            show.legend = FALSE,
            size = 3,
            parse = TRUE) +
  geom_text(data = lin_mods,
            aes(x = x,
                y = y3,
                label = pval_text),
            show.legend = FALSE,
            size = 3) +
  xlab("Median daily water table (m)") + 
  ylab("Understory species count or Shannon's H")
ggsave(veghyd,
       filename = "veg_models/div_wt.tiff",
       device = "tiff",
       dpi = 300)
