# 
# Author: Jake Diamond
# Purpose: Analysis of surface models, distributions, semivariograms
# Date: February 4, 2019
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(MASS)
library(gstat)
library(tidyverse)

# Run analysis for different detrended data types
for(i in 1:3){
  # Get all filenames
  if(i == 1){
    print("Doing analysis for quadtratic trended data")
    filenames <- paste("Lidar/detrended_quad/", 
                       list.files("Lidar/detrended_quad"), 
                       sep = "")
    trend <- "quad"
  } else if(i ==2){
    print("Doing analysis for linear detrended data")
    filenames <- paste("Lidar/detrended/", 
                       list.files("Lidar/detrended"), 
                       sep = "")
    trend <- "linear"
  } else{
    print("Doing analysis for non-detrended data")
    filenames <- paste("Lidar/Rasters_density_filtered/", 
                       list.files("Lidar/Rasters_density_filtered"), 
                       sep = "")
    trend <- "no_trend"
  }
  # rm(p_hum_all, z_l_all)
  # Bimodal analysis --------------------------------------------------------
  for (j in 1:length(filenames)){
    
    # Get site from file name
    s <- filenames[j] %>%
      strsplit(., "[/]") %>%
      unlist()
    s <- s[[3]] %>%
      strsplit(., "[_]") %>%
      unlist()
    s <- s[[1]]

    # Load raster file
    r <- raster(filenames[j])
    
    # Convert raster to points
    xyz <- rasterToPoints(r)
    xyz <- as.data.frame(xyz)
    colnames(xyz)[3] <- "z"
    
    # Get rid of data outliers (like super high points that are probably branches)
    xyz <- dplyr::filter(xyz, z < quantile(z, 0.99)
                         , z > quantile(z, 0.001)
    )
    # Normalize elevations by location at the well
    well <- read.csv("relative_elevations_all_v6.csv",
                     stringsAsFactors = FALSE) %>%
      dplyr::filter(site == s,
                    point == "well") %>%
      dplyr::select(z, z_mod_lin, z_mod_quad)
      
    well <- ifelse(trend == "no_trend",
                   well$z,
                   ifelse(trend == "linear",
                          well$z - well$z_mod_lin,
                          well$z - well$z_mod_quad))
    xyz$z <- xyz$z - well
    
    # Remove data for memory
    rm(r)
    
    # First want to get rid of edge effects, so reduce the
    # window size
    xyzSub <- xyz %>%
      dplyr::filter(between(x, 
                            quantile(x, 0.05),
                            quantile(x, 0.95)),
                    between(y, 
                            quantile(y, 0.05),
                            quantile(y, 0.95)))

    # Semivariogram analysis --------------------------------------------------
  #   print(paste("Doing semivariogram analysis for site", s))
    # Turn data in to a spatial dataframe for analysis
    xy <- xyzSub[, c(1, 2)]
    spdf <- SpatialPointsDataFrame(coords = xy,
                                   data = xyzSub,
                                   proj4string =
                                     CRS("+proj=utm +zone=15 +datum=WGS84"))

    # Subsample data for faster computation
    spSample <- spdf[sample(1:length(spdf), 20000), ]

    # # basic variogram
    vario <- variogram(z~1,
                       data = spSample,
                       width = 1)
  # 
  #   # directional variogram
  #   vario_dir <- variogram(z~1, data = spSample,
  #                          alpha = c(0, 45, 90, 135))
  #   
  #   # Plot directional variogram
  #   jpeg(paste(s, trend, "directional_semivariogram.jpg",
  #              sep = "_"))
  #   print(plot(vario_dir,
  #        xlab = "Distance (m)",
  #        ylab = "Semivariance"))
  #   dev.off()
  #   
  #   # Estimate best variogram model fit
    v_fit <- fit.variogram(vario, 
                           vgm(model = c("Mat", "Exc", "Exp",
                                         "Nug", "Ste")), 
                               fit.kappa = TRUE)
    v_fit$site <- s
    # model info
    # v_fit
    # plot(vario, v_fit)
    # mo <- v_fit[1, 1]
    zsill <- var(spSample$z)
    # nug <- v_fit[1, 2]
    # ran <- v_fit[2, 3]
    # 
    # # Plot data
    # semi <- ggplot(data = vario, aes(x = dist, y = gamma)) +
  #     geom_point(size = 1.8) + theme_bw() + 
  #     theme(legend.position = "none",
  #           panel.grid.major.x = element_blank(),
  #           panel.grid.minor.x = element_blank(),
  #           panel.grid.minor.y = element_blank(),
  #           panel.grid.major.y = element_blank()) +
  #     geom_hline(yintercept =  zsill, 
  #                size = 1.4,
  #                linetype = "dashed") + 
  #     annotate("text", x = 2.5, y = zsill + 0.002, 
  #              label = paste("Sill =", 
  #                            round(zsill, 3))) +
  #     ylab("Semivariance") +
  #     xlab("Distance (m)")
  #   semi
  #   ggsave(plot = semi, 
  #          filename = paste(s, 
  #                           trend,
  #                           "semivariogram.tiff",
  #                           sep = "_"),
  #          device = "tiff",
  #          width = 4, height = 3, 
  #          units = "in")
  #   
    # Collect all semivariogram data in one place
    semi_site <- data.frame(dist = vario$dist,
                            gam = vario$gamma,
                            site = s,
                            sill = zsill)
    if(j == 1){
      semi_results <- semi_site
      semi_fits <- v_fit
    } else{
      semi_results <- bind_rows(semi_results, semi_site)
      semi_fits <- bind_rows(semi_fits, v_fit)
    }
  }

  print("Writing semivariogram analysis to .csv")
  write.csv(semi_results, file = paste0("semi_1m_results_",
                                        trend, ".csv"))
  write.csv(semi_fits, file = paste0("semi_1m_fit_results_",
                                     trend, ".csv"))
}
