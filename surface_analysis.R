# 
# Author: Jake Diamond
# Purpose: Analysis of surface models, distributions, semivariograms
# Date: December 1, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)

# Get raster files from point cloud models
r <- raster("Lidar/Clouds/T1_1cm.asc")

# Read in xy coordinates
df_xy <- data.frame(x = c(1, 2, 3), 
                    y = c(1, 2, 3)
                    )
# Extract raster values from xy coordinates
xy_r_values <- extract(r, 
                       SpatialPoints(xy), 
                       sp = TRUE
                       )@data 

# 2) you want to analyze the distributions in what way? If you just want the raw numbers to throw in a model of the size distribution just grab the Z values (raw) or the scalar field (detrended) of .asc files. You would read it like this, I believe:
  
  dist<-data.frame(Z = read.table("my_model.asc", sep = ",")[,3], Zd =  read.table("my_model.asc", sep = ",")[,7])
  #from here you can run any stats you'd like to analyze the distribution