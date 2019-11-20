# 
# Author: Jake Diamond
# Purpose: To analyze black ash understory vegetation community matrix with NMDS
# Date: February 26, 2019
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(MASS)
library(tidyverse)
library(broom)
library(vegan)
library(ggrepel)
# ggplot theme
theme_set(theme_bw() + 
            theme(panel.grid = element_blank()))
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
  group_by(site, hu.ho, species) %>%
  dplyr::summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0) %>%
  rename_all(funs(make.names(.))) %>%
  dplyr::select(-starts_with("Forb"),
                -starts_with( "Unknown"),
                -"lightly.spiky") %>%
  arrange(site, hu.ho) %>%
  ungroup()

# Get groups of data
groups <- data.frame(site = df_nmds$site,
                     hu_ho = df_nmds$hu.ho)

# rarefy
set.seed(42)
min_samp <- round(min(rowSums(df_nmds[, -(1:2)])))
df_nmds_rare <- as.data.frame(round(rrarefy(df_nmds[, -(1:2)], min_samp)))

# Dissimilarity matrix
nmds_distr <- as.matrix((vegdist(df_nmds_rare, "bray")))
nmds_dist <- as.matrix((vegdist(df_nmds[,-(1:2)], "bray")))
#perform NMDS
NMDSr <- metaMDS(nmds_distr)
stressplot(NMDSr)
plot(NMDSr)
NMDSr
NMDS <- metaMDS(nmds_dist)
# Testing
adonis_huho <- adonis(nmds_distr ~ hu_ho, groups)
adonis_huho
adonis_huho_s <- adonis(nmds_distr ~ hu_ho, 
                        strata = groups$site, groups)
adonis_huho_s
adonis_site <- adonis(nmds_distr ~ site, groups)
adonis_site

#build a data frame with NMDS coordinates and metadata
MDS1 <- NMDSr$points[,1]
MDS2 <- NMDSr$points[,2]
NMDSp <- data.frame(MDS1 = MDS1, 
                  MDS2 = MDS2, 
                  site = groups$site, 
                  huho = groups$hu_ho)

# Plot
NMDS_p <- ggplot(NMDSp, 
       aes(x = MDS1, y = MDS2, color = huho)) +
  geom_text(aes(label = site),
            show.legend = FALSE) +
  scale_color_manual(name = "Microsite",
                      labels = c("Hollow", "Hummock"),
                      values = c("grey50", "black")) +
  stat_ellipse() +
  coord_fixed() + 
  theme(legend.position = c(0.14, 0.87),
               legend.background = element_rect(fill = "grey90",
                                                color = "black")) + 
  ylab("NMDS1") +
  xlab("NMDS2")
NMDS_p
ggsave(NMDS_p,
       filename = "veg_models/NMDS_vegetation_bw.tiff",
       device = "tiff",
       dpi = 300)



