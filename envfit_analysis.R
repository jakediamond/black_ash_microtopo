# 
# Author: Jake Diamond
# Purpose: To analyze understory community data with vegan pkg
# Date: February 22, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(MASS)
library(tidyverse)
library(broom)
library(vegan)
library(ggrepel)

# Load veg data
df <- readRDS("veg_data.RDS")
df$species <- as.character(df$species)
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"
# Quick rename
df[df$species == "Aster lateriflorus (L) Britton", "species"] <- 
  "Aster lateriflorus"

# Get community matrix with sites-hu.ho-pt as index
df_nmds <- df %>%
  group_by(site, hu.ho, point, species) %>%
  dplyr::summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0) %>%
  rename_all(funs(make.names(.))) %>%
  dplyr::select(-starts_with("Forb"),
                -starts_with( "Unknown"),
                -"lightly.spiky") %>%
  arrange(site, point) %>%
  ungroup() %>%
  dplyr::select(-site, -point, -hu.ho)

# Load elevation data
elev <- readRDS("elevations_wt") %>%
  arrange(site, point) %>%
  dplyr::select(site, point, z_relh_mean) %>%
  group_by(site, point) %>%
  summarize_all(funs(mean), na.rm = T)

# Load soil chem data
chem <- readRDS("soil_chem_wide") %>%
  dplyr::select(site, point, ca, so4, mg, no3, po4, cl, 
                C, CN, N, z, z_well, hu.ho) %>%
  mutate(z = z - z_well) %>%
  dplyr::select(-z_well)

# Scale function
scale_this <- function(x){
  (x /  mean(x, na.rm=TRUE))
}

# Relative concentration
chema <- chem %>%
  group_by(site) %>%
  mutate_at(funs(scale_this), 
            .vars = vars(-point, -site, -hu.ho))
  

# Combine to get all environmental data
ba_env <- left_join(elev, chem, by = c("site", "point")) %>%
  dplyr::select(-point)
# Determine which rows have no env data
noenv <- which(!complete.cases(ba_env))

# Remove data from the species and env matrix (clipped data)
df_nmds_sub <- df_nmds[-noenv, ]
ba_env <- ba_env[-noenv,]

# Get rid of rare species (found in <1% total samples)
# First save which species these are for methods
species_rare <- df_nmds_sub[, 
                        colSums(df_nmds != 0) <= 
                          0.01 * nrow(df_nmds)]
df_nmds_sub <- df_nmds_sub[,
                       colSums(df_nmds != 0) > 
                         0.01 * nrow(df_nmds)]
# Get groups
groups <- data.frame(site = ba_env$site, huho = ba_env$hu.ho)
# NMDS
set.seed(42)
nmds <- metaMDS(df_nmds_sub, trace = F,
                k = 4,
                trymax = 100)
nmds
plot(nmds, type = "t")
stressplot(nmds)
# Get scores for plotting in ggplot2
scrs <- as.data.frame(scores(nmds, display = "sites"))
scrs.sp <- as.data.frame(scores(nmds, display = "species"))
scrs <- cbind(scrs,
              # species = rownames(scrs)
              groups
              )
scrs.sp <- cbind(scrs.sp,
                 species = rownames(scrs.sp))
# Environmental fitting
ef <- envfit(nmds, ba_env[, -c(1, 2, 13)], permu = 999,
             na.rm = TRUE)
ef
env.scrs <- as.data.frame(scores(ef, display = "vectors"))
env.scrs <- cbind(env.scrs, env = rownames(env.scrs),
                  pval = ef$vectors$pvals)
# Rename environmental variables for plotting
env.scrs$env <- c("SO[4]^{`2-`}",
                  "Mg^{`2+`}",
                  "NO[3]^{`-`}",
                  "PO[4]^{`3-`}",
                  "C",
                  "CN",
                  "N",
                  "z")

# Plotting
nmdsp <- ggplot(scrs) +
  geom_point(aes(x = NMDS1, y = NMDS2,
                           color = site,
                 shape = huho),
             size = 2) +
  coord_fixed() +
  geom_segment(data = env.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey60") +
  geom_text(data = env.scrs,
            aes(x = NMDS1, y = NMDS2, label = env),
            size = 4,
            parse = TRUE) +
  scale_color_viridis_d(name = "Site") +
  scale_shape_manual(values = c(1, 16),
                     breaks = c("hollow", "hummock"),
                     labels = c("Hollow", "Hummock"),
                     name = "Microsite") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.12, 0.4),
        legend.background = element_rect(fill = "grey90",
                                         color = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 2))
nmdsp
ggsave(nmdsp,
       filename = "veg_models/nmds_environment_v2.tiff",
       dpi = 300,
       device = "tiff")
# Dissimilarities and environment
# Analyse vegetation–environment relationships without
# ordination, or in full space. adonis 
# implements a multivariate analysis of variances using distance matrices
# We define beta diversity as the  slope of species-area curve, 
# or the exponent z of the Arrhenius model where the number of species 
# S is dependent on the size X of the study area. For pairwise 
# comparison of sites the slope z can be found from the
# number of species shared between two sites (a) and the 
# number of species
# unique to each sites (b and c). It is commonly regarded that z ≈ 0.3
# implies random sampling variability, and only higher values mean real
# systematic differences. The Arrhenius z can be directly found with function betadiver that also provided many other indices of pairwise beta
# diversity
betad <- betadiver(df_nmds_sub, "z")
plot(betad)
# Function adonis can use formula interface, and the dependent data can
# be either dissimilarities or data frame, and in the latter case adonis uses
# vegdist to find the dissimilarities.
adonis(betad ~ site*hu.ho, ba_env, perm = 300)

# Function adonis studied the differences in the group means, but function betadisper studies the differences in group homogeneities. Function
# adonis was analogous to multivariate analysis of variance, and betadisper is analogous to Levene’s test of the equality of variances
mod <- with(ba_env, betadisper(betad, hu.ho))
mod
plot(mod)
boxplot(mod)
anova(mod)
permutest(mod)
TukeyHSD(mod)


# Mantel test
# Mantel test compares two sets of dissimilarities
# In this example we study how well the lichen pastures (varespec) correspond to the environment. We have already used vector fitting after ordination. However, the ordination and environment may be non-linearly
# related, and we try now with function mantel. We first perform a pca of
# environmental variables, and then compute dissimilarities for first principal components. We use standard R function prcomp, but princomp or
# rda will work as well. Function scores in vegan will work with all these
# methods. The following uses the same standardizations for community
# dissimilarities as previously used in metaMDS
pc <- prcomp(ba_env[,-c(1,2,13)], scale = TRUE)
pc <- scores(pc, display = "sites", choices = 1:4)
edis <- vegdist(pc, method = "euclid")
vare.dis <- vegdist(wisconsin(sqrt(df_nmds_sub)))
mantel(vare.dis, edis)
plot(vare.dis, edis)
# Everything is O.K. if the relationship is more or less monotonous, or even
# linear and positive. In spatial models we even may observe a hump which
# indicates spatial aggregation.