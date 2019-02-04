# 
# Author: Jake Diamond
# Purpose: Detrend surface models to remove broad slopes for more consistent comparison across points
# Date: December 1, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(tidyverse)

# Get all filenames of .tif raster files for xyz data of sites
filenames <- paste("Lidar/Rasters_density_filtered/", 
                   list.files("Lidar/Rasters_density_filtered"), 
                   sep = "")

# Load raw elevation data of measurement points
elev <- read.csv("Lidar/Point_validation/valpts_r.csv",
                 stringsAsFactors = FALSE)
elev$site <- toupper(elev$site)
elev[elev$site == "B1", "site"] <- "L1"
elev[elev$site == "B3", "site"] <- "L2"
elev[elev$site == "B6", "site"] <- "L3"

# Load raw elevation data of delineated hummocks
hums <- read.csv("hummock_stats_clean.csv",
                 stringsAsFactors = FALSE)
hums$X <- NULL
hums <- hums %>%
  dplyr::rename(x = xmean,
                y = ymean,
                z = zmean)
hums[hums$site == "B1", "site"] <- "L1"
hums[hums$site == "B3", "site"] <- "L2"
hums[hums$site == "B6", "site"] <- "L3"

# Combine data
elev <- bind_rows(elev, hums)

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
  fit_linear <- lm(z ~ x + y, data = xyz)
  
  # Fit a quadratic detrend
  fit_quad <- lm(z ~ poly(x, y, 
                          degree = 2),
                 data = xyz)
  
  # Get residuals. These are the corrections to raw elevations
  resid_linear <- fit_linear$residuals
  resid_quad <- fit_quad$residuals
  
  # Recreate the dataframe with coordinates and residuals
  xyz$resid_lin <- resid_linear
  xyz$resid_quad <- resid_quad
  # Remove data for less memory
  # rm(resid.linear)
  rm(resid_quad, resid_linear)
  
  # convert back to rasterlayer
  # Rasterize based on original raster extent, 
  # possibly used for spatial analysis
  # Only need to do this the first time running
  # r_detrend_lin <- rasterize(xyz[, 1:2], 
  #                            r, 
  #                            xyz$resid_lin, 
  #                            fun = mean)
  # names(r_detrend_lin) <- "z"
  # print(paste("Writing detrended rasters for site", s))
  # # Write the detrended raster to disc
  # writeRaster(r_detrend_lin,
  #             filename = paste0(s, 
  #                               "_all_1cm_linear_detrend.tif"),
  #             format = "GTiff",
  #             datatype = "FLT8S",
  #             options = c("COMPRESS=NONE",
  #                         'TFW = YES'),
  #             overwrite = TRUE
  #             )
  # # Same for quadratic detrend data
  # r_detrend_quad <- rasterize(xyz[, 1:2], r, 
  #                             xyz$resid_quad, fun = mean)
  # # Write the detrended raster to disc
  # writeRaster(r_detrend_quad,
  #             filename = paste0(s, "_all_1cm_quad_detrend.tif"),
  #             format = "GTiff",
  #             datatype = "FLT8S",
  #             options = c("COMPRESS=NONE",
  #                         'TFW = YES'),
  #             overwrite = TRUE
  #             )
  # Remove data for less memory
  rm(xyz, r_detrend_lin, r_detrend_quad, r)
  # Calculate detrended data
  detrended <- dplyr::filter(elev,
                             site == s) %>%
    dplyr::select(x, y) %>%
    mutate(z_mod_lin = predict(fit_quad, newdata = .),
           z_mod_quad = predict(fit_linear, newdata = .)) %>%
    inner_join(elev, by = c("x", "y"))
  rm(fit_linear, fit_quad)
  print(paste("Writing detrended elevation points for site", s))
  if(i == 1){
    result <- detrended
  } else {
    result <- rbind(result, detrended)
  }
}
# write.csv(result, "valpts_detrend_r.csv")
write.csv(result, "clean_hummocks_and_valpts_detrend.csv")
