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
data_path_h <- "Lidar/Hummock_delineation/Cleaned"
files_h <- dir(data_path_h)
# Get data path for rasters
data_path_r <- "Lidar/Rasters/"
files_r <- dir(data_path_r, pattern = "cm.tif")

# Read in all hummock clean hummock files and merge
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
                    ~ read_csv(file.path(data_path_h, .)))) %>%
  dplyr::select(-filename) %>%
  rename(hum = data) %>%
  unnest() %>%
  rename(x = `//X`, y = Y, id = treeID) %>%
  dplyr::select(-Z,-Classification, -`Original cloud index`)

# rasterxyz function
# rxyzfun <- function(data){
#   data %>% 
#     rasterToPoints() %>%
#     as.data.frame() %>%
#     rename(x = 1, y = 2)
# }
# Get all rasters to xy data to get convex hulls for holes in polygons
# hums_df <- hums %>%
#   mutate(xyz = map(hum, rxyzfun)) %>%
#   dplyr::select(-hum) %>%
#   unnest() %>%
#   mutate(id = coalesce(B1_humID, B3_humID, D1_humID, 
#                        D3_humID, D4_humID, T1_humID)) %>%
#   dplyr::select(site, x, y, id)
# Convex hulls
ch <- hums %>% 
  nest(-c(site, id)) %>%
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
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_D3 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "D3") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_D4 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "D4") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_L1 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "L1") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_L2 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "L2") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}
pch_T1 <- ch %>%
  mutate(id = as.character(id)) %>%
  dplyr::filter(site == "T1") %>%
  dplyr::select(-site) %>%
  nest(-id) %>%
  mutate(Poly = map(data, ~select(., x, y) %>% Polygon()),
         polys = map2(Poly, id, ~Polygons(list(.x),.y))) %>%
         {SpatialPolygons(.$polys)}

# Merge those
hums <- tibble(site = sites,
               hum = c(pch_D1,pch_D3, pch_D4,
                       pch_L1, pch_L2, pch_T1))

# Read in all surface model raster files and merge
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
humstats <- read_rds("hummock_stats_clean_detrended_z") %>%
  mutate(zmaxn = ifelse(site %in% c("L1", "L2", "L3"),
                        zmeann,
                        zmaxn),
         zmeann = ifelse(site %in% c("L1", "L2", "L3"),
                         z20n,
                         zmeann))

# Combine with hummock ids
df <- left_join(match, humstats, by = c("site", "id"))

# Normalize elevations by location at the well
df <- read.csv("relative_elevations_all_v6.csv",
               stringsAsFactors = FALSE) %>%
  dplyr::filter(
    point == "well") %>%
  dplyr::select(site, z) %>%
  rename(zwell = z) %>%
  right_join(df) %>%
  mutate(zrel = z - zwell,
         zuse = ifelse(is.na(zmax),
                       z,
                       zmax),
         zuse = zuse - zwell,
         zuse = ifelse(site == "L1",
                       z,
                       zmax))

# number of matches per site


# Quick plot
ggplot(df, aes(x = zmaxn, y = dbh, color = site)) + 
  geom_point() + geom_smooth(method = "glm") +
  facet_wrap(~site, scales = "free_x") 
write_rds(df, "tree_dbh_hum_match")
df2 <- read_rds("tree_dbh_hum_match")
