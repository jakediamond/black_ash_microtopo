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
library(viridis)

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
  rm(p_hum_all, z_l_all)
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
    
    # Make sure names align
    s_hu <- ifelse(s == "L1",
                   "B1",
                   ifelse(s == "L2",
                          "B3",
                          ifelse(s == "L3",
                                 "B6",
                                 s)))
    
    # Get hummock data for plotting
    hu <- read.csv(paste0("Lidar/Hummock_delineation/Cleaned/",
                          s_hu,
                          ".txt")) %>%
      mutate(Classification = 1) %>%
      rename(X = X..X, hum = Classification, id = treeID) %>%
      left_join(read.csv("hummock_stats_clean.csv",
                         stringsAsFactors = FALSE) %>%
                  dplyr::filter(site == s_hu) %>%
                  dplyr::select(-X),
                by = c("id")) %>%
      dplyr::filter(area > 0.1,
                    zmaxn > 0.15) %>%
      dplyr::select(X, Y, Z, hum, id) %>%
      mutate(X = round(X, 2),
             Y = round(Y, 2))
    # Get data into long format and subsample for fast proc.
    z_l <- xyz %>%
      mutate(X = round(x, 2),
             Y = round(y, 2)) %>%
      left_join(hu, by = c("X", "Y")) %>%
      mutate(hum = ifelse(!is.na(hum), 1,
                          ifelse(z < mean(z),
                                 0,
                                 NA))) %>%
      dplyr::filter(!is.na(hum)) %>%
      sample_n(10000) %>%
      dplyr::select(z, hum)

    dif <- mean(z_l[z_l$hum == 1, "z"]) - 
      mean(z_l[z_l$hum == 0, "z"])
    dif_sd <- sqrt(sd(z_l[z_l$hum == 1, "z"])^2 + 
                     sd(z_l[z_l$hum == 0, "z"])^2)
    
    # Plot data colored by hummock or hollow
    p_hum <- ggplot() +
      geom_density(data = z_l,
                     aes(x = z,
                         y=(..scaled..) * n,
                         fill = as.factor(hum)),
                   alpha = 0.3) + 
      scale_fill_viridis(breaks = c(0, 1),
                        labels = c("Hollow",
                                   "Hummock"),
                        discrete = TRUE,
                        option = "cividis") +
      theme_bw() +
      # geom_text(x = 0.15,
      #           y = 0.65,
      #           aes(label = paste("list(Delta*bar(z) ==", 
      #                             round(dif, digits=2), ")")),
      #           show.legend = FALSE,
      #           parse = TRUE) +
      theme(legend.position = c(0.75, 0.8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank()
      ) +
      ylab("Density") +
      xlab("Relative Elevation (m)")
    p_hum
    ggsave(plot = p_hum,
           filename = paste(s, 
                            trend,
                            "histogram.tiff",
                            sep = "_"),
           device = "tiff",
           width = 4, height = 3,
           units = "in")
      
    
    # Plot bimodal if D is less than -6 (rule of thumb)
    # And also if the difference in means is greater than 10 cm
    # if(D > 20 & D_m > 0.1){
    #   mean1 <- bimdl$parameters$mean[1]
    #   mean2 <- bimdl$parameters$mean[2]
    #   sigma1 <- sqrt(bimdl$parameters$variance$sigmasq[1])
    #   sigma2 <- ifelse(is.na(sqrt(bimdl$parameters$variance$sigmasq[2])),
    #                    sigma1,
    #                    sqrt(bimdl$parameters$variance$sigmasq[2]))
    #   lam1 <- bimdl$parameters$pro[1]
    #   lam2 <- bimdl$parameters$pro[2]
    #   # Function to plot both normal distributions
    #   plot_mix_comps <- function(x, mu, sigma, lam) {
    #     lam * dnorm(x, mu, sigma)
    #   }
    #   # Plotting bimodal
    #   p_b <- ggplot() +
    #     geom_histogram(data = xyz,
    #                    aes(x = z,
    #                        ..density..,
    #                        fill =..x..),
    #                    bins = 100) +
    #     stat_function(fun = plot_mix_comps,
    #                   aes(colour = "1"),
    #                   args = list(mean1,
    #                               sigma1,
    #                               lam1),
    #                   lwd = 1.2) +
    #     stat_function(fun = plot_mix_comps,
    #                   aes(colour = "2"),
    #                   args = list(mean2,
    #                               sigma2,
    #                               lam2),
    #                   lwd = 1.2) +
    #     geom_segment(aes(x = mean1,
    #                      xend = mean1,
    #                      y = 0,
    #                      yend = max(density(xyz$z)$y)),
    #                  colour = "blue",
    #                  linetype = "dashed",
    #                  size = 1.2
    #     ) +
    #     geom_segment(aes(x = mean2,
    #                      xend = mean2,
    #                      y = 0,
    #                      yend = max(density(xyz$z)$y)),
    #                  colour = "green",
    #                  linetype = "dashed",
    #                  size = 1.2
    #     ) +
    #     geom_text(aes(x = mean1 + 0.08,
    #                   y = max(density(xyz$z)$y) + 0.05,
    #                   label = round(mean1, 2)),
    #               color = "blue"
    #     ) +
    #     geom_text(aes(x = mean2 + 0.08,
    #                   y = max(density(xyz$z)$y) + 0.05,
    #                   label = round(mean2, 2)),
    #               color = "green"
    #     ) +
    #     scale_colour_manual("Component",
    #                         values = c("1" = "blue",
    #                                    "2" = "green")) +
    #     scale_fill_gradientn(colours = c("blue",
    #                                      "green",
    #                                      "yellow",
    #                                      "red")) +
    #     theme_bw() +
    #     theme(legend.position = "none",
    #           axis.title.x = element_text(face = "bold",
    #                                       vjust = 0.6),
    #           axis.title.y = element_text(face = "bold",
    #                                       vjust = 0.6),
    #           # axis.text.x = element_text(size = 24),
    #           # axis.text.y = element_text(size = 24),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank()
    #     ) +
    #     ylab("Density") +
    #     xlab("Relative Elevation (m)")
    #   
    #   ggsave(plot = p_b, filename = paste(s, 
    #                                       trend,
    #                                       "Bimodal.tiff",
    #                                       sep = "_"),
    #          device = "tiff",
    #          width = 4, height = 3,
    #          units = "in")
    # } else {
    #   mean1 <- bimdl$parameters$mean[1]
    #   sigma1 <- sqrt(bimdl$parameters$variance$sigmasq[1])
    #   # Function to plot normal distributions
    #   plot_comps <- function(x, mu, sigma, lam) {
    #     dnorm(x, mu, sigma)
    #   }
    #   # Plotting unimodal
    #   p_u <- ggplot(data = xyz, aes(x = z)) +
    #     # ,fill = ..x..)) +
    #     geom_histogram(aes(z,
    #                        ..density..),
    #                    bins = 100) +
    #     stat_function(fun = plot_comps,
    #                   args = list(mean1,
    #                               sigma1),
    #                   lwd = 1.2,
    #                   color = "black") +
    #     # scale_fill_gradientn(colours = c("blue",
    #     #                                  "green",
    #     #                                  "yellow",
    #     #                                  "red")) +
    #     theme_bw() +
    #     theme(legend.position = "none",
    #           axis.title.x = element_text(face = "bold",
    #                                       vjust = 0.6),
    #           axis.title.y = element_text(face = "bold",
    #                                       vjust = 0.6),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank()
    #     ) +
    #     ylab("Density") +
    #     xlab("Relative Elevation (m)")
    #   
    #   ggsave(plot = p_u, filename = paste(s, 
    #                                       trend,
    #                                       "_Unimodal.tiff",
    #                                       sep = "_"),
    #          device = "tiff",
    #          width = 4, height = 3,
    #          units = "in")
    # }
    # Combine all cluster data for all sites
    mcl_site <- mcl_tidy %>%
      dplyr::rename(G = component) %>%
      full_join(mcl_bic, by = c("site","G"))
    
    z_l$site <- s
    
    if(j == 1){
      mcl_results <- mcl_site
      z_l_all <- z_l
    } else{
      mcl_results <- bind_rows(mcl_results, mcl_site)
      z_l_all <- bind_rows(z_l_all, z_l)
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
               sep = "_"))
    print(plot(vario_dir,
         xlab = "Distance (m)",
         ylab = "Semivariance"))
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
  
  print("Plotting semivariogram analysis")
  # Remove a bunch of data
  rm(bimdl, dmcl, hu, mcl, mcl_bic, mcl_results, mcl_site,
     mcl_tidy, p_hum, semi, semi_fits, semi_site, spdf,
     spSample, trimdl, unimdl, v_fit, vario, vario_dir, xy, xyz,
     z_l, xyzSample)
  # Write the z_l_data to disc
  write_rds(z_l_all,
            path = paste0("hummock_elevations_",
                          trend))
  
  # Add a column for plotting about site type
  semi_results <- semi_results %>%
    mutate(type = str_sub(site, start = 1L, end = 1L))
  # Plot all data on one plot for comparison
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
    # geom_hline(yintercept =  zsill, 
    #            size = 1.2,
    #            linetype = "dashed") + 
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
  # Remove data
  rm(semi_all, semi_results)
  
  # Plot data colored by hummock or hollow
  p_hum_all <- ggplot() +
    geom_density(data = z_l_all,
                 aes(x = z,
                     fill = as.factor(hum)),
                 alpha = 0.3) + 
    scale_fill_viridis(breaks = c(0, 1),
                       labels = c("Hollow",
                                  "Hummock"),
                       discrete = TRUE,
                       option = "cividis") +
    theme_bw() +
    # geom_text(x = 0.15,
    #           y = 0.65,
    #           aes(label = paste("list(Delta*bar(z) ==", 
    #                             round(dif, digits=2), ")")),
    #           show.legend = FALSE,
    #           parse = TRUE) +
    theme(legend.position = c(0.85, 0.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()
    ) +
    facet_wrap(~site) +
    ylab("Density") +
    xlab("Relative Elevation (m)")
  
  ggsave(plot = p_hum_all,
         filename = paste("all_sites", 
                          trend,
                          "histogram.tiff",
                          sep = "_"),
         device = "tiff",
         width = 6, height = 4,
         units = "in")
}

# Plot all data based on best detrend
# L3 no detrend
# L2 linear
# L1 no detrend
# D1 quad (linear??)
# D2 linear
# D3 no detrend
# D4 quad
# T1 quad
# T2 no detrend
# T3 no detrend
# No detrend sites
nd <- c("L3", "L1", "D3", "T2", "T3")
l <- c("L2", "D2")
q <- c("D1", "D4", "T1")


# Get all data based on best detrend
df <- read_rds("hummock_elevations_no_trend") %>%
  dplyr::filter(site %in% nd) %>%
  bind_rows(read_rds("hummock_elevations_linear") %>%
              dplyr::filter(site %in% l)) %>%
  bind_rows(read_rds("hummock_elevations_quad") %>%
              dplyr::filter(site %in% q))

# Summarize data
df_sum <- read.csv("Lidar/hummock_stats_ext6.csv") %>%
  dplyr::select(id, site_area) %>%
  mutate(site = str_sub(id, 1, 2)) %>%
  dplyr::select(site, site_area) %>%
  distinct() %>%
  right_join(read.csv("hummock_stats_clean.csv",
                   stringsAsFactors = FALSE)) %>%
  dplyr::select(-X) %>%
  dplyr::filter(area > 0.1,
                zmaxn > 0.15) %>%
  group_by(site) %>%
  summarize(no = n(),
            aratio = sum(area) / mean(site_area),
            zavg = mean(zmaxn),
            zsd = sd(zmaxn)) %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site))))

df$type <- str_sub(df$site, 1, 1)
df$num <- str_sub(df$site, 2, 2)
df_sum$type <- str_sub(df_sum$site, 1, 1)
df_sum$num <- str_sub(df_sum$site, 2, 2)

# Get x and y coordinates for the plot
df_sum$x <- c()

# Subsample data for faster processing
dfSam <- dplyr::sample_n(df, size = 10000)

library(ggpubr)
p_hum_best_d <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "D"),
               aes(x = z,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  theme_bw() +
  # geom_text(data = df_sum,
  #           x = x,
  #           y = y,
  #           aes(label = paste("list(Delta*bar(z) ==",
  #                             round(zavg, digits=2), 
  #                             "%+-%", round(zsd, digits=2), ")")),
  #           show.legend = FALSE,
  #           parse = TRUE,
  #           size = 2.5) +
  # geom_text(data = df_sum,
  #           x = x,
  #           y = y,
  #           aes(label = paste("list(n ==",
  #                             no, ")")),
  #           show.legend = FALSE,
  #           parse = TRUE,
  #           size = 2.5) +
  # geom_text(data = df_sum,
  #           x = x,
  #           y = y,
  #           aes(label = paste("list(A[ratio] ==",
  #                             round(aratio, digits=2), ")")),
  #           show.legend = FALSE,
  #           parse = TRUE,
  #           size = 2.5) +
  theme(legend.justification = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_da <- p_hum_best_d + theme(legend.position = "none")

p_hum_best_l <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "L"),
               aes(x = z,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  theme_bw() +
  # geom_text(data = df_sum,
  #           x = x,
  #           y = y,
  #           aes(label = paste("list(Delta*bar(z) ==",
  #                             round(zavg, digits=2), 
  #                             "%+-%", round(zsd, digits=2), ")")),
  #           show.legend = FALSE,
  #           parse = TRUE,
  #           size = 2.5) +
  # geom_text(data = df_sum,
  #           x = x,
#           y = y,
#           aes(label = paste("list(n ==",
#                             no, ")")),
#           show.legend = FALSE,
#           parse = TRUE,
#           size = 2.5) +
# geom_text(data = df_sum,
#           x = x,
#           y = y,
#           aes(label = paste("list(A[ratio] ==",
#                             round(aratio, digits=2), ")")),
#           show.legend = FALSE,
#           parse = TRUE,
#           size = 2.5) +
theme(legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      axis.text.y = element_blank()
) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")
p_hum_best_l

p_hum_best_t <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "T"),
               aes(x = z,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  theme_bw() +
  # geom_text(data = df_sum,
  #           x = x,
  #           y = y,
  #           aes(label = paste("list(Delta*bar(z) ==",
  #                             round(zavg, digits=2), 
  #                             "%+-%", round(zsd, digits=2), ")")),
  #           show.legend = FALSE,
  #           parse = TRUE,
  #           size = 2.5) +
  # geom_text(data = df_sum,
  #           x = x,
#           y = y,
#           aes(label = paste("list(n ==",
#                             no, ")")),
#           show.legend = FALSE,
#           parse = TRUE,
#           size = 2.5) +
# geom_text(data = df_sum,
#           x = x,
#           y = y,
#           aes(label = paste("list(A[ratio] ==",
#                             round(aratio, digits=2), ")")),
#           show.legend = FALSE,
#           parse = TRUE,
#           size = 2.5) +
theme(legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      axis.text.y = element_blank()
) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")
p_hum_best_t

legend <- cowplot::get_legend(p_hum_best_d)
library(cowplot)
p_hum_best2 <- ggdraw() +
  draw_plot(p_hum_best_da + rremove("x.text") + rremove("x.title"), 
            x = 0, y = 0.7, width = 1, height = 0.3) +
  draw_plot(p_hum_best_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.4, width = 0.76, height = 0.3) +
  draw_plot(p_hum_best_t,
            x = 0, y = 0, width = 0.76, height = 0.4) + 
  draw_grob(legend,
            0.82, 0, .3/3.3, 0.5)
p_hum_best2
ggsave(plot = p_hum_best2,
       filename = "all_sites_densities_v2.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")
