# 
# Author: Jake Diamond
# Purpose: Find lowest points of detrended rasters
# Date: December 7, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(tidyverse)

# Get all filenames of .tif raster files for xyz data of sites
filenames <- paste("Lidar/Rasters/detrended/", 
                   list.files("Lidar/Rasters/detrended"), 
                   sep = "")

# Load raw elevation data of measurement points
elev <- read.csv("valpts_detrend_r.csv",
                 stringsAsFactors = FALSE)


# Loop through every site and find minimum elevation of detrended surface 
# append detrended data to an elevation datatable associated with measured points.
for (i in 1:length(filenames)){
  # Get site from file name
  s <- filenames[i] %>%
    strsplit(., "[/]") %>%
    unlist()
  s <- s[[4]] %>%
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
  # Minimum point
  minpt <- data.frame(site = s,
                      minpt = min(xyz$z))
  if(i == 1){
    result <- minpt
  } else {
    result <- rbind(result, minpt)
  }
}

elev <- left_join(elev, result)
elev$X <- NULL
elev$site <- as.factor(elev$site)
write_rds(elev, "elevation_data")
