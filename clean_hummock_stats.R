# 
# Author: Jake Diamond
# Purpose: Hummock statistics for cleaned data
# Date: February 1, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(rgeos)
library(sp)

# Get filenames of all cleaned hummock files
data_path <- "Lidar/Hummock_delineation/Cleaned"
files <- dir(data_path, pattern = "*.txt") 

# Read in all files and merge
df <- tibble(filename = files,
                 site = str_sub(filename, 1, 2)) %>% 
  mutate(data = map(filename,
                    ~ read_csv(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  rename(x = `//X`, y = Y, z = Z, id = treeID) %>%
  dplyr::select(1:5)

# Hummock statistics function
humstats_fun <- function (data, res){
  colnames(data)[1:3]<-c("X","Y","Z")
  humxy <- data.frame(x = data[, 1], y = data[, 2])
  # Convert to polygon for areas with holes
  hum_p <- SpatialPolygons(list(
    Polygons(
      list(
        Polygon(
          humxy[chull(humxy[, 1], humxy[, 2]), ]
        )
      ), ID = 1
    )
  ),
  proj4string = 
    CRS("+proj=utm +zone=15 +datum=WGS84")
  )
  # Tools to calculate on polygons, accounts for 
  # tree trunk holes but messes up on weird hummock shapes
  cent <- gCentroid(hum_p)
  centroid_x <- cent@coords[1]
  centroid_y <- cent@coords[2]
  area_poly <- gArea(hum_p)
  perim_poly <- gLength(hum_p)
  PA_poly <- perim_poly / area_poly
  return(data.frame(zmin = min(data$Z),
                    zmax = max(data$Z),
                    zmean = mean(data$Z),
                    xmin = min(data$X),
                    xmax = max(data$X),
                    xmean = mean(data$X),
                    ymin = min(data$Y),
                    ymax = max(data$Y),
                    ymean = mean(data$Y),
                    zmaxn = max(data$Z - min(data$Z)),
                    z20 = quantile(data$Z, 0.2),
                    z80 = quantile(data$Z, 0.8),
                    z20n = quantile(data$Z - min(data$Z), 0.2),
                    z80n = quantile(data$Z - min(data$Z), 0.8),
                    zmeann = mean(data$Z - min(data$Z)),
                    area = res^2*nrow(data),
                    vol = sum((data$Z - min(data$Z))*(res^2)),
                    centx = centroid_x,
                    centy = centroid_y,
                    area_poly = area_poly,
                    perim_poly = perim_poly,
                    PA_poly = PA_poly)
         )
} 

# analysis of hummock stats
df_stat <- df %>%
  group_by(site, id) %>%
  nest() %>%
  mutate(stats = map(data, humstats_fun, res = 0.01)) %>%
  unnest(stats, .drop = TRUE)

# Write data to disc
write.csv(df_stat, "hummock_stats_clean.csv")

ggplot(data =filter(df_stat, area > 0.1),
       aes(x = z20n)) +geom_histogram()+facet_wrap(~site)
