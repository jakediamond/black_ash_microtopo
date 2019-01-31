# 
# Author: Jake Diamond
# Purpose: Do density filter on point clouds and save as rasters
# Date: January 26, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(tidyverse)
library(lidR)
library(data.table)

# Get all filenames
# Note: T1 is riddled with hummocks, not looking bimodal
# D2, L3 maybe needs to be detrended (quad)
filenames <- paste("Lidar/1_cm_update/", 
                   list.files("Lidar/1_cm_update"), 
                   sep = "")
# Bimodal analysis --------------------------------------------------------
for (i in 1:length(filenames)){
  # Get site from file name
  s <- filenames[i] %>%
    strsplit(., "[/]") %>%
    unlist()
  s <- s[[3]] %>%
    strsplit(., "[_]") %>%
    unlist()
  s <- s[[1]]
  
  # Load point cloud
  cloud <- fread(filenames[i],
                select = c(1:3), 
                skip = 0, 
                header = FALSE, 
                sep = " ")
  colnames(cloud) <- c("X", "Y", "Z")
  las <- LAS(cbind(cloud[, 1:3]))
  rm(cloud)
  writeLAS(las, 
           paste(
                 s, 
                 "_cloud_raw.las", 
                 sep = "")
           )
  
  # Resample the data based on return densities > 2000 pts/m2
  print(paste("Reading and filtering LAS surface model from site", s))
  las <- readLAS(paste0(s, "_cloud_raw.las"))
  density <- grid_density(las, res = 1)
  rm(las)
  r <- raster(paste0("Lidar/Rasters/", 
                     ifelse(s == "B1",
                            "L1",
                            ifelse(s == "B3",
                                   "L2",
                                   ifelse(s == "B6",
                                          "L3",
                                          s))),
                     "_1cm.tif"))
  density <- as.raster(density)
  density <- resample(density, r)
  r[density < 2000] <- NA
  writeRaster(r,
              filename = paste0(ifelse(s == "B1",
                                       "L1",
                                       ifelse(s == "B3",
                                              "L2",
                                              ifelse(s == "B6",
                                                     "L3",
                                                     s))), 
                                "_density.tif"),
              format = "GTiff",
              datatype = "FLT8S",
              options = c("COMPRESS=NONE",
                          'TFW = YES'),
              overwrite = TRUE)
}
