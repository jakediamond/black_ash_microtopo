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
library(mclust)
library(MASS)
library(gstat)
library(tidyverse)
library(broom)

# Note: T1 is riddled with hummocks, not looking bimodal
# D2, L3 maybe needs to be detrended (quad)
# D1, L2 (maybe), T2, T3 (maybe) needs to be detrended (linear)
# D3, D4, L3 does not need to be detrended
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
  # Bimodal analysis --------------------------------------------------------
  for (j in 1:length(filenames)){
    # rm(p_b, p_u)
    # Get site from file name
    s <- filenames[j] %>%
      strsplit(., "[/]") %>%
      unlist()
    s <- s[[3]] %>%
      strsplit(., "[_]") %>%
      unlist()
    s <- s[[1]]
    print(paste("Doing cluster analysis for site", s))
    # Load raster file
    r <- raster(filenames[j])
    # plot(r)
    # Convert raster to points
    xyz <- rasterToPoints(r)
    xyz <- as.data.frame(xyz)
    colnames(xyz)[3] <- "z"
    
    # Get rid of data outliers (like super high points that are probably branches)
    xyz <- dplyr::filter(xyz, z < quantile(z, 0.99)
                         , z > quantile(z, 0.001)
    )
    # Normalize elevations by the lowest elevations
    xyz$z <- xyz$z - quantile(xyz$z, 0.001)
    
    # Remove data for memory
    rm(r)
    
    # Subsample the data for less computing power...same results
    xyzSample <- dplyr::sample_n(xyz, size = 10000)
    
    # Get an esimate of the number of clusters that defines the sample
    mcl <- Mclust(xyzSample$z)
    dmcl <- densityMclust(xyzSample$z)
    mclustBIC(xyzSample$z)
    
    # Tidy the data for export
    mcl_tidy <- tidy(mcl) %>%
      dplyr::filter(row_number()==n()) %>%
      dplyr::select(5:ncol(.)) %>%
      gather(key = "component",
             value = "mean") %>%
      mutate(component = as.numeric(str_sub(component, 
                                            start = -1L))) %>%
      right_join(tidy(mcl) %>%
                   dplyr::select(1:4)) %>%
      mutate(site = s)
    
    # Compare bimodal, unimodal, and trimodal
    unimdl <- Mclust(xyzSample$z,
                     G = 1)
    bimdl <- Mclust(xyzSample$z,
                    G = 2)
    trimdl <- Mclust(xyzSample$z,
                     G = 3)
    
    # Get data frame of BIC values
    mcl_bic <- glance(mcl) %>%
      bind_rows(glance(unimdl),
                glance(bimdl),
                glance(trimdl)) %>%
      mutate(site = s) %>%
      arrange(G) %>%
      distinct(G, model, .keep_all = TRUE)
    
    print(paste("Plotting cluster analysis for site", s))
    # Plot best fit
    jpeg(paste(s, trend, "best_fit_density_model.jpg",
               sep = "_"), 
         width = 800, height = 800)
    plot(dmcl,
         what = "density",
         data = xyzSample$z,
         breaks = 30,
         xlab="Relative elevation (m)",
         main = paste(s, "density"))
    dev.off()
    
    # calculate BIC difference between bi- and unimodal
    D <- bimdl$bic - unimdl$bic
    
    # Calculate mean difference in bimodal dist
    D_m <- bimdl$parameters$mean[2] - bimdl$parameters$mean[1]
    
    # Plot bimodal if D is less than -6 (rule of thumb)
    # And also if the difference in means is greater than 10 cm
    if(D > 20 & D_m > 0.1){
      mean1 <- bimdl$parameters$mean[1]
      mean2 <- bimdl$parameters$mean[2]
      sigma1 <- sqrt(bimdl$parameters$variance$sigmasq[1])
      sigma2 <- ifelse(is.na(sqrt(bimdl$parameters$variance$sigmasq[2])),
                       sigma1,
                       sqrt(bimdl$parameters$variance$sigmasq[2]))
      lam1 <- bimdl$parameters$pro[1]
      lam2 <- bimdl$parameters$pro[2]
      # Function to plot both normal distributions
      plot_mix_comps <- function(x, mu, sigma, lam) {
        lam * dnorm(x, mu, sigma)
      }
      # Plotting bimodal
      p_b <- ggplot() +
        geom_histogram(data = xyz,
                       aes(x = z,
                           ..density..,
                           fill =..x..),
                       bins = 100) +
        stat_function(fun = plot_mix_comps,
                      aes(colour = "1"),
                      args = list(mean1,
                                  sigma1,
                                  lam1),
                      lwd = 1.2) +
        stat_function(fun = plot_mix_comps,
                      aes(colour = "2"),
                      args = list(mean2,
                                  sigma2,
                                  lam2),
                      lwd = 1.2) +
        geom_segment(aes(x = mean1,
                         xend = mean1,
                         y = 0,
                         yend = max(density(xyz$z)$y)),
                     colour = "blue",
                     linetype = "dashed",
                     size = 1.2
        ) +
        geom_segment(aes(x = mean2,
                         xend = mean2,
                         y = 0,
                         yend = max(density(xyz$z)$y)),
                     colour = "green",
                     linetype = "dashed",
                     size = 1.2
        ) +
        geom_text(aes(x = mean1 + 0.08,
                      y = max(density(xyz$z)$y) + 0.05,
                      label = round(mean1, 2)),
                  color = "blue"
        ) +
        geom_text(aes(x = mean2 + 0.08,
                      y = max(density(xyz$z)$y) + 0.05,
                      label = round(mean2, 2)),
                  color = "green"
        ) +
        scale_colour_manual("Component",
                            values = c("1" = "blue",
                                       "2" = "green")) +
        scale_fill_gradientn(colours = c("blue",
                                         "green",
                                         "yellow",
                                         "red")) +
        theme_bw() +
        theme(legend.position = "none",
              axis.title.x = element_text(face = "bold",
                                          vjust = 0.6),
              axis.title.y = element_text(face = "bold",
                                          vjust = 0.6),
              # axis.text.x = element_text(size = 24),
              # axis.text.y = element_text(size = 24),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        ) +
        ylab("Density") +
        xlab("Relative Elevation (m)")
      
      ggsave(plot = p_b, filename = paste(s, 
                                          trend,
                                          "Bimodal.tiff",
                                          sep = "_"),
             device = "tiff",
             width = 4, height = 3,
             units = "in")
    } else {
      mean1 <- bimdl$parameters$mean[1]
      sigma1 <- sqrt(bimdl$parameters$variance$sigmasq[1])
      # Function to plot normal distributions
      plot_comps <- function(x, mu, sigma, lam) {
        dnorm(x, mu, sigma)
      }
      # Plotting unimodal
      p_u <- ggplot(data = xyz, aes(x = z)) +
        # ,fill = ..x..)) +
        geom_histogram(aes(z,
                           ..density..),
                       bins = 100) +
        stat_function(fun = plot_comps,
                      args = list(mean1,
                                  sigma1),
                      lwd = 1.2,
                      color = "black") +
        # scale_fill_gradientn(colours = c("blue",
        #                                  "green",
        #                                  "yellow",
        #                                  "red")) +
        theme_bw() +
        theme(legend.position = "none",
              axis.title.x = element_text(face = "bold",
                                          vjust = 0.6),
              axis.title.y = element_text(face = "bold",
                                          vjust = 0.6),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        ) +
        ylab("Density") +
        xlab("Relative Elevation (m)")
      
      ggsave(plot = p_u, filename = paste(s, 
                                          trend,
                                          "_Unimodal.tiff",
                                          sep = "_"),
             device = "tiff",
             width = 4, height = 3,
             units = "in")
    }
    # Combine all cluster data for all sites
    mcl_site <- mcl_tidy %>%
      dplyr::rename(G = component) %>%
      full_join(mcl_bic, by = c("site","G"))
    
    if(j == 1){
      mcl_results <- mcl_site
    } else{
      mcl_results <- bind_rows(mcl_results, mcl_site)
    }
    
    # Semivariogram analysis --------------------------------------------------
    print(paste("Doing semivariogram analysis for site", s))
    # Turn data in to a spatial dataframe for analysis
    xy <- xyz[, c(1, 2)]
    spdf <- SpatialPointsDataFrame(coords = xy,
                                   data = xyz,
                                   proj4string =
                                     CRS("+proj=utm +zone=15 +datum=WGS84"))
    
    # Subsample data for faster computation
    spSample <- spdf[sample(1:length(spdf), 10000), ]
    
    # basic variogram
    vario <- variogram(z~1, 
                       data = spSample, 
                       width = 0.5)

    # directional variogram
    vario_dir <- variogram(z~1, data = spSample,
                           alpha = c(0, 45, 90, 135))
    
    # Plot directional variogram
    jpeg(paste(s, trend, "directional_semivariogram.jpg",
               sep = "_"), 
         width = 800, height = 800)
    plot(vario_dir,
         xlab = "Distance (m)",
         ylab = "Semivariance")
    dev.off()
    
    # Estimate best variogram model fit
    v_fit <- fit.variogram(vario, vgm(c("Exp", 
                                        "Gau", 
                                        "Sph")))
    v_fit$site <- s
    # Sill of data
    zsill <- var(spSample$z)
    
    # Plot data
    semi <- ggplot(data = vario, aes(x = dist, y = gamma)) +
      geom_point(size = 1.8) + theme_bw() + 
      theme(legend.position = "none",
            axis.title.x = element_text(face = "bold", 
                                        vjust = 0.6), 
            axis.title.y = element_text(face = "bold", 
                                        vjust = 0.6),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank()) +
      geom_hline(yintercept =  zsill, 
                 size = 1.4,
                 linetype = "dashed") + 
      # geom_vline(xintercept = locmin,
      #            size = 1,
      #            linetype = "dashed",
      #            color = "dark grey") +
      annotate("text", x = 2.5, y = zsill + 0.002, 
               label = paste("Sill =", 
                             round(zsill, 3))) +
      # scale_x_continuous(breaks = seq(0, 15, 2)) +
      # scale_y_continuous(breaks = seq(0, 0.03, 0.005)) +
      ylab("Semivariance") +
      xlab("Distance (m)")
    
    ggsave(plot = semi, 
           filename = paste(s, 
                            trend,
                            "semivariogram.tiff",
                            sep = "_"),
           device = "tiff",
           width = 4, height = 3, 
           units = "in")
    
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
  print("Writing cluster analysis to .csv")
  write.csv(mcl_results, file = paste0("mcl_results_", 
                                       trend, ".csv"))
  
  print("Writing semivariogram analysis to .csv")
  write.csv(semi_results, file = paste0("semi_results_", 
                                        trend, ".csv"))
  write.csv(semi_fits, file = paste0("semi_fit_results_", 
                                     trend, ".csv"))
  
  print("Plotting semivagiogram analysis")
  # Add a column for plotting about site type
  semi_results <- semi_results %>%
    mutate(type = str_sub(site, start = 1L, end = 1L))
  # Plot all data on one plto for comparison
  semi_all <- ggplot(data = semi_results, 
                     aes(x = dist, 
                         y = gam,
                         colour = site,
                         shape = type)) +
    scale_shape_manual(breaks = c("D", "L", "T"),
                       values = c(16, 1, 0),
                       guide = FALSE) + 
    scale_colour_viridis_d(name = "Site") + 
    # guide_legend(nrow = 2) +
    geom_point(size = 1.8) + theme_bw() + 
    theme(legend.position = "bottom",
          axis.title.x = element_text(face = "bold", 
                                      vjust = 0.6), 
          axis.title.y = element_text(face = "bold", 
                                      vjust = 0.6),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    geom_hline(yintercept =  zsill, 
               size = 1.2,
               linetype = "dashed") + 
    # geom_vline(xintercept = locmin,
    #            size = 1,
    #            linetype = "dashed",
    #            color = "dark grey") +
    # annotate("text", x = 2.5, y = zsill + 0.002, 
    #          label = paste("Sill =", 
    #                        round(zsill, 3))) +
    # scale_x_continuous(breaks = seq(0, 15, 2)) +
    # scale_y_continuous(breaks = seq(0, 0.03, 0.005)) +
    ylab("Semivariance") +
    xlab("Distance (m)")
  
  ggsave(plot = semi_all, 
         filename = paste(trend,
                          "all_semivariograms.tiff",
                          sep = "_"),
         device = "tiff",
         width = 4, height = 3, 
         units = "in")
}
