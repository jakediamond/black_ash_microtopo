# 
# Author: Jake Diamond
# Purpose: To analyze soil chemistry data
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/Soils/Soil Extraction Chemistry")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(PerformanceAnalytics)
library(tidyverse)
library(broom)
library(lubridate)
library(viridis)
library(corrplot)
library(plot3D)

# Load chemistry data
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

# Linear model for elevation against solutes
lm_model <- function(data) {
  z = data$z - data$z_well
  lm(conc_rel ~ z, data = data)
}

# Site types and numbers
df$type <- str_sub(df$site, 1, 1)
df$num <- str_sub(df$site, 2, 2)

# Linear analyses of confining depth vs elevation
lin_mods_site <- df %>%
  group_by(site, type, num, solute) %>%
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

lin_mods <- df %>%
  group_by(solute) %>%
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
# Correlation plot
x <- df %>%
  ungroup() %>%
  spread(solute, conc_rel) %>%
  dplyr::select(C:so4, z) %>%
  # mutate(z_rel = z - mean,
  #        zd_rel = zd - mean) %>%
  ungroup() %>%
  dplyr::select(C:so4, z) %>%
  group_by(z) %>%
  summarize_all(funs(mean), na.rm = T) %>%
  dplyr::select(-z) %>%
  as.matrix()

# Save this correlation matrix for later to compare with veg data
# write_rds(x, "soil_chem_matrix")
# Plot of soil chem correlations
png("corrmat_chem_global_rel.png",
    width = 1000,
    height = 1000)
chart.Correlation(x, histogram=TRUE, pch=19, method = "spearman")
dev.off()
# Strip panel names for bar plot
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

# Plot
chem_p <- ggplot(data = df,
                                      # !(solute == "Cl^{`-`}" &
                                      #     conc > 30)
                                      # , !(solute == "NO[3]^{`-`}" &
                                      #     conc > 1)
                                      # , !(site %in% c("L1",
                                      #               "L2",
                                      #               "L3"))
                                      # , site %in% c("D3", "T1",
                                      #               "T2", "D4")
                                      # , log == 0
                                      # ),
                 aes(x = z - z_well,
                     y = conc_rel
                     # , color = site
                     )
                 ) +
  geom_point(aes(shape = hu.ho)) +
  geom_smooth(method = "lm", se = FALSE)+
  # geom_smooth(se = FALSE) +
  scale_shape_manual(name = "",
                     values = c(1, 16)) +
  scale_colour_viridis(discrete = TRUE) +
  ylab(expression("Soil Extraction Concentration (mg "*L^{-1}*")")
       ) +
  xlab("Relative elevation (m)") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~solute, scales = "free",
             labeller = label_parsed)
chem_p


# 3D plots of all soil data -----------------------------------------------
sites <- unique(df_data_w$site)
solutes <- c("ca", "cl", "mg", "no3", 
             "po4", "so4", "C", "N", "CN")
for(i in 1:length(sites)){
  sit = sites[i]
  df_3d <- dplyr::filter(df_data_w, site == sit)
  for(j in 1:length(solutes)){
    solute = solutes[j]
    color_pts = dplyr::pull(df_3d, solute)
    tit = paste(sit, solute)
    png(filename = paste(
      sit,
      solute,
      "quad.png", 
      sep = "_"), 
      width = 800, 
      height = 800)
    text3D(df_3d$x, 
              df_3d$y, 
              df_3d$z_rel.d, 
              bty = "b2", 
              pch = 19,
              main = tit,
              labels = df_3d$point,
              colvar = color_pts, 
              cex = 2, 
              ticktype = "detailed"
              , theta = 45,
              phi = 30
              # ,
              # colkey = list(at = c(2, 3, 4), side = 1, 
              # addlines = TRUE, 
              # length = 0.5, 
              # width = 0.5
              # ,labels = c("setosa", 
              # "versicolor", 
              # "virginica")
              # ) 
    )
    dev.off()
    
  }
}