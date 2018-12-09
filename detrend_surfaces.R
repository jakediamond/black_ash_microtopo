# 
# Author: Jake Diamond
# Purpose: Detrend surface models to 
# Date: December 1, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(tidyverse)

# Get all filenames of .tif raster files for xyz data of sites
filenames <- paste("Lidar/Rasters/", 
                   list.files("Lidar/Rasters"), 
                   sep = "")

# Load raw elevation data of measurement points
elev <- read.csv("Lidar/Point_validation/valpts_r.csv",
                 stringsAsFactors = FALSE)
elev$site <- toupper(elev$site)
elev[elev$site == "B1", "site"] <- "L1"
elev[elev$site == "B3", "site"] <- "L2"
elev[elev$site == "B6", "site"] <- "L3"

# Loop through every site and detrend, save results as raster, and also 
# append detrended data to an elevation datatable associated with measured points.
for (i in 1:length(filenames)){
  # Get site from file name
  s <- filenames[i] %>%
    strsplit(., "[/]") %>%
    unlist()
  s <- s[[3]] %>%
    strsplit(., "[_]") %>%
    unlist()
  s <- s[[1]]
  # Load raster file
  r <- raster(filenames[i])
  # Convert raster to points data
  xyz <- rasterToPoints(r)
  # Convert points data to a dataframe
  xyz <- as.data.frame(xyz)
  colnames(xyz)[3] <- "z"
  
  # Fit a linear detrend
  fit.linear <- lm(z ~ x + y, data = xyz)
  
  # Get residuals. These are the corrections to raw elevations
  resid.linear <- fit.linear$residuals
  
  # Recreate the dataframe with coordinates and residuals
  xyz$resid <- resid.linear
  # Remove data for less memory
  rm(resid.linear)
  
  # convert back to rasterlayer
  # Rasterize based on original raster extent, possibly used for spatial analysis
  # Only need to do this the first time running
  # r.detrend <- rasterize(xyz[, 1:2], r, xyz[, 4], fun = mean)
  # writeRaster(r.detrend,
  #             paste0(s, "_1cm_detrend.tif"),
  #             options = c('TFW = YES')
  # )
  # Remove data for less memory
  rm(xyz, r.detrend, r)
  # Calculate detrended data
  detrended <- dplyr::filter(elev,
                             site == s) %>%
    select(x, y) %>%
    mutate(z_mod = predict(fit.linear, newdata = .)) %>%
    inner_join(elev, by = c("x", "y"))
  if(i == 1){
    result <- detrended
  } else {
    result <- rbind(result, detrended)
  }
}

write.csv(result, "valpts_detrend_r.csv")
