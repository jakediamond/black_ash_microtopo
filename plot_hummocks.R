# 
# Purpose: To plot delineated hummocks for every black ash site
# Author: Jake Diamond
# Date: February 1, 2019
# 
# Set working directory
setwd("C:/Users/diamo/Dropbox/Projects/EAB/")

# Load libraries
library(ggplot2)
library(ggspatial)
library(raster)
library(sf)
library(scales)
library(viridis)  # better colors for everyone

# Load data
r <- raster("delineation_code/output/D1_humID.tif")
r_site <- raster("Data/Lidar/Rasters_density_filtered/D1_density.tif")

# Convert data to dataframe
r_spdf <- as(r, "SpatialPixelsDataFrame")
r_df <- as.data.frame(r_spdf)
colnames(r_df) <- c("id", "x", "y")

rs_spdf <- as(r_site, "SpatialPixelsDataFrame")
rs_df <- as.data.frame(rs_spdf)
colnames(rs_df) <- c("x", "y", "z")

# ggplot() + 
#   layer_spatial(r_spdf, aes(fill = D1_humID))

# Plot data
ggplot() +  
  # geom_tile(data = r_df, aes(x = x, y = y, fill = id), alpha = 0.8) + 
  geom_tile(data = rs_df, aes(x = x, y = y, fill = NA), alpha = 0.4) + 
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"))
