# 
# Author: Jake Diamond
# Purpose: To analyze soil chemistry data
# Date: September 3, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data/Soils/Soil Extraction Chemistry")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data/Soils")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)
library(viridis)
library(corrplot)
library(plot3D)

# Load data
df <- read.csv("Soil Extraction Chemistry/meta_for_r.csv",
               stringsAsFactors = FALSE)
df_chem <- read.csv("Soil Extraction Chemistry/soil_chem_r.csv",
                    stringsAsFactors = FALSE)
df_cn <- read.csv("CN_r.csv")
# elev <- read.csv("C:/Users/diamo/Dropbox/Projects/EAB/Data/relative_elevations.csv")
elev <- read.csv("C:/Users/diamo/Dropbox/Projects/EAB/Data/relative_elevations_quad.csv")
elev$X <- NULL
df <- read.csv("depth_to_confining.csv",
                    stringsAsFactors = FALSE) %>%
  fill(plot) %>%
  right_join(df)

# Make sure that depth to confining layer is numeric
# and account for data >120 cm
df$depth <- ifelse(df$depth == ">120",
                   "150",
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

# Correlation plot
x <- df_data_w %>%
  dplyr::filter(sample_id != 196) %>%
  mutate(z_rel_mean = z_rel - mean,
         z_rel.d_mean = z_rel.d - mean,
         z_rel_conf = z_rel - depth / 100,
         z_rel.d_conf = z_rel.d - depth / 100,
         z_conf = z - depth / 100) %>%
  dplyr::select(-sample_id,
         -(site:intermediate),
         -point_id,
         -(mean:quant_20),
         -(z_well:z_well.d)) %>%
  as.matrix() %>%
  na.omit()
corrplot(cor(x), method="number")

# Strip panel names for plot
df_data_l$solute <- factor(df_data_l$solute,
                           levels = c("ca",
                                      "cl",
                                      "mg",
                                      "no3",
                                      "po4",
                                      "so4",
                                      "C",
                                      "N",
                                      "CN"))
levels(df_data_l$solute) <- c("Ca^{`2+`}",
                         "Cl^{`-`}",
                         "Mg^{`2+`}",
                         "NO[3]^{`-`}",
                         "PO[4]^{`3-`}",
                         "SO[4]^{`2-`}",
                         "'%'*C",
                         "'%'*N",
                         "C:N")

# Plot
chem_p <- ggplot(data = dplyr::filter(df_data_l,
                                      !(solute == "Cl^{`-`}" &
                                          conc > 30)
                                      , !(solute == "NO[3]^{`-`}" &
                                          conc > 1)
                                      # , !(site %in% c("L1",
                                      #               "L2",
                                      #               "L3"))
                                      , site %in% c("D3", "T1",
                                                    "T2", "D4")
                                      , log == 0
                                      ),
                 aes(x = z_rel.d - quant_20,
                     y = conc
                     , color = site
                     )
                 ) +
  geom_point(aes(shape = hu.ho)) +
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


