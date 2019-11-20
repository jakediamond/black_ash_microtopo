# 
# Author: Jake Diamond
# Purpose: To perform indicator species analysis on black ash sites
# Date: February 26, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(indicspecies)

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
  summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0) %>%
  rename_all(funs(make.names(.))) %>%
  dplyr::select(-starts_with("Forb"),
                -starts_with( "Unknown"),
                -"lightly.spiky") %>%
  arrange(site, point) %>%
  ungroup()

# Get groups
groups <- df_nmds[,1:2]
df_nmds <- df_nmds[, -(1:3)]

# Indicator species on site
indval_s <- multipatt(df_nmds, groups$site,
                    control = how(nperm=999))
# Summary of significant components 
summary(indval_s, indvalcomp = TRUE)

# Indicator species analysis on hummock hollow
indval <- multipatt(df_nmds, groups$hu.ho,
                    control = how(nperm=999))
# Summary of significant components 
summary(indval, indvalcomp = TRUE)
# Summary of all components
summary(indval, alpha = 1)
# presence absence ("r.g" accoutns for some groups having more sites than others)
pa <- as.data.frame(ifelse(df_nmds>0,1,0))
phi <- multipatt(pa, groups$hu.ho, func = "r.g",
                 control = how(nperm=999))
summary(phi)

# We can reduce the number of combinations to consider by selecting the 
# candidate species to combine. 
# In this example, we choose those species
# whose frequency within Group 2 is larger than 40%
B=strassoc(df_nmds, cluster=groups$hu.ho ,func="B")
sel=which(B[,1]>0.1)
sel
summary(B)
# indicators should in most cases be run with the option 
# requesting for bootstrap confidence intervals 
# (using option nboot)
sc= indicators(X=df_nmds[,sel], 
               cluster=groups$hu.ho, group="hollow", 
               verbose=TRUE,
               At=0.5, Bt=0.2,
               nboot = 200)
# print the results of indicators for the most useful indicators
print(sc)

sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
