# 
# Purpose: Get base elevations of trees in black ash wetlands
# Author: Jake Diamond
# Date: February 18, 2019
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(tidyverse)
library(broom)

# Load xy data of dbh by site
df <- read.csv("Lidar/Trees/All_DBH_xymatch.csv")
sites <- as.character(unique(df$site))
# Get data path for hummock ids
data_path_h <- "E:/Dropbox/Dropbox/Projects/EAB/delineation_code/output"
files_h <- dir(data_path_h, pattern = "ID.tif")
# Get data path for rasters
data_path_r <- "Lidar/Rasters/"
files_r <- dir(data_path_r, pattern = "cm.tif")

# Read in all hummock raster files and merge
hums <- tibble(filename = files_h,
               site = str_sub(filename, 1, 2)) %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site)))) %>% 
  dplyr::filter(site %in% sites) %>%
  mutate(data = map(filename,
                    ~ raster(file.path(data_path_h, .)))) %>%
  dplyr::select(-filename) %>%
  rename(hum = data)

# rasterxyz function
rxyzfun <- function(data){
  data %>% 
    rasterToPoints() %>%
    as.data.frame() %>%
    rename(x = 1, y = 2)
}
# Get all rasters to xy data to get convex hulls for holes in polygons
hums_df <- hums %>%
  mutate(xyz = map(hum, rxyzfun)) %>%
  dplyr::select(-hum) %>%
  unnest() %>%
  mutate(id = coalesce(B1_humID, B3_humID, D1_humID, 
                       D3_humID, D4_humID, T1_humID)) %>%
  dplyr::select(site, x, y, id)
# Convex hulls
ch <- hums_df %>% 
  group_by(site, id) %>%
  nest() %>%
  mutate(
    hull = map(data, ~ with(.x, chull(x, y))),
    out = map2(data, hull, ~ .x[.y,,drop=FALSE])
  ) %>%
  dplyr::select(-data, -hull) %>%
  unnest()

# Polygons of convex hulls for spatial match
pch_D1 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "D1") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_D3 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "D3") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_D4 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "D4") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_L1 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "L1") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_L2 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "L2") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_T1 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "T1") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~dplyr::select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}

# Merge those
hums <- tibble(site = sites,
               hum = c(pch_D1, pch_D3, pch_D4,
                       pch_L1, pch_L2, pch_T1))

# Read in all hummock raster files and merge
rs <- tibble(filename = files_r,
               site = str_sub(filename, 1, 2)) %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site)))) %>% 
  dplyr::filter(site %in% sites) %>%
  mutate(data = map(filename,
                    ~ raster(file.path(data_path_r, .)))
  ) %>%
  dplyr::select(-filename)

# Get all points to be spatial points
pts <- df %>%
  group_by(site) %>%
  nest() %>%
  mutate(xy = map(data, dplyr::select, 
                  x, y),
         xy = map(xy, SpatialPoints)) %>%
  rename(dbh = data)

# Get spatial matches of points onto hummock rasters
match <- pts %>%
  left_join(hums) %>%
  left_join(rs) %>%
  group_by(site) %>%
  mutate(matchums = map2(xy, hum, sp::over),
         matchraster = map2(data, xy, raster::extract)) %>%
  dplyr::select(-xy, -data, -hum) %>%
  unnest() %>%
  dplyr::select(-z) %>%
  rename(z = matchraster, id = matchums)

# Load hummock stats
humstats <- read.csv("Lidar/hummock_stats_ext3.csv",
                     stringsAsFactors = FALSE)
humstats <- humstats %>%
  dplyr::select(-X) %>%
  separate(id, into = c("site", "id"), "_")
humstats$id <- as.numeric(humstats$id)
humstats[humstats$site == "B1", "site"] <- "L1"
humstats[humstats$site == "B3", "site"] <- "L2"
humstats[humstats$site == "B6", "site"] <- "L3"
# Combine with hummock ids
df <- left_join(match, humstats, by = c("site", "id"))

# Normalize elevations by location at the well
df2 <- read.csv("relative_elevations_all_v6.csv",
               stringsAsFactors = FALSE) %>%
  dplyr::filter(
    point == "well") %>%
  dplyr::select(site, z) %>%
  rename(zwell = z) %>%
  right_join(df) %>%
  mutate(
         zuse = ifelse(is.na(zmin_raw),
                       z,
                       zmin_raw),
         zuse = zuse- zwell) 

# write to disc
write_rds(df2, "dbh_hummock_match_minz")
