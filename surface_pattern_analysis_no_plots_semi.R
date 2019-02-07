# 
# Author: Jake Diamond
# Purpose: Analysis of surface models, distributions, semivariograms
# Date: February 4, 2019
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
library(cowplot)
library(ggpubr)

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
    print(paste("Doing cluster analysis for site", s))
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
                            quantile(x, 0.1),
                            quantile(x, 0.9)),
                    between(y, 
                            quantile(y, 0.1),
                            quantile(y, 0.9)))
    # Subsample the data for less computing power...same results
    # xyzSample <- dplyr::sample_n(xyzSub, size = 10000)
    # 
    # # Get an esimate of the number of clusters that defines the sample
    # mcl <- Mclust(xyzSample$z)
    # dmcl <- densityMclust(xyzSample$z)
    # mclustBIC(xyzSample$z)
    # 
    # # Tidy the data for export
    # mcl_tidy <- tidy(mcl) %>%
    #   dplyr::filter(row_number()==n()) %>%
    #   dplyr::select(5:ncol(.)) %>%
    #   gather(key = "component",
    #          value = "mean") %>%
    #   mutate(component = as.numeric(str_sub(component, 
    #                                         start = -1L))) %>%
    #   right_join(tidy(mcl) %>%
    #                dplyr::select(1:4)) %>%
    #   mutate(site = s)
    # 
    # # Compare bimodal, unimodal, and trimodal
    # unimdl <- Mclust(xyzSample$z,
    #                  G = 1)
    # bimdl <- Mclust(xyzSample$z,
    #                 G = 2)
    # trimdl <- Mclust(xyzSample$z,
    #                  G = 3)
    # 
    # # Get data frame of BIC values
    # mcl_bic <- glance(mcl) %>%
    #   bind_rows(glance(unimdl),
    #             glance(bimdl),
    #             glance(trimdl)) %>%
    #   mutate(site = s) %>%
    #   arrange(G) %>%
    #   distinct(G, model, .keep_all = TRUE)
    # 
    # print(paste("Plotting cluster analysis for site", s))
    # # Plot best fit
    # jpeg(paste(s, trend, "best_fit_density_model.jpg",
    #            sep = "_"), 
    #      width = 800, height = 800)
    # plot(dmcl,
    #      what = "density",
    #      data = xyzSample$z,
    #      breaks = 30,
    #      xlab="Relative elevation (m)",
    #      main = paste(s, "density"))
    # dev.off()
    # 
    # # calculate BIC difference between bi- and unimodal
    # D <- bimdl$bic - unimdl$bic
    # 
    # # Calculate mean difference in bimodal dist
    # D_m <- bimdl$parameters$mean[2] - bimdl$parameters$mean[1]
    
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
    # Join the hummock metadata with the xyz data
    z_l <- xyzSub %>%
      mutate(X = round(x, 2),
             Y = round(y, 2)) %>%
      left_join(hu, by = c("X", "Y")) %>%
      mutate(hum = ifelse(!is.na(hum), 1, 0),
             hum2 = ifelse(hum == 1, 1,
                           ifelse(z < mean(z),
                                 0,
                                 NA))) %>%
      dplyr::select(-X, -Y, -Z, -x, -y)

    # Subsample data for faster processing
    z_lSam <- z_l %>%
      sample_n(10000)

    # Write stuff for smaller files
    write_rds(z_l,
              path = paste0("RDS_files/",
                            s, "_hummock_elevations_",
                            trend))
    write_rds(z_lSam,
              path = paste0("RDS_files/Samples/",
                            s, "_hummock_elevations_sample",
                            trend))
    
    # Remove data
    rm(hu, dmcl, mcl, bimdl, trimdl, unimdl)
    
    # Plot data colored by hummock or hollow
    # p_hum <- ggplot(data = z_lSam) +
    #   geom_density(aes(x = z,
    #                    y = (..scaled..) * n,
    #                    fill = as.factor(hum)),
    #                alpha = 0.3) + 
    #   scale_fill_viridis(breaks = c(0, 1),
    #                      labels = c("Hollow",
    #                                 "Hummock"),
    #                     discrete = TRUE,
    #                     option = "cividis") +
    #   theme_bw() +
    #   # geom_text(y = 4200,
    #   #           aes(x = max(z) - (max(z)- min(z)) * 0.25,
    #   #               label = paste("list(Delta*bar(z) ==",
    #   #                             round(dif, digits=2),
    #   #                             "%+-%", 
    #   #                             round(dif_sd, digits=2), 
    #   #                             ")")),
    #   #           show.legend = FALSE,
    #   #           parse = TRUE) +
    #   theme(legend.position = c(0.75, 0.8),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         legend.title = element_blank()
    #   ) +
    #   ylab("Count") +
    #   xlab("Relative Elevation (m)")
    # p_hum
    # ggsave(plot = p_hum,
    #        filename = paste(s, 
    #                         trend,
    #                         "histogram_sample.tiff",
    #                         sep = "_"),
    #        device = "tiff",
    #        width = 4, height = 3,
    #        units = "in")
     
    # Combine all cluster data for all sites
    # mcl_site <- mcl_tidy %>%
    #   dplyr::rename(G = component) %>%
    #   full_join(mcl_bic, by = c("site","G"))
    # 
    # z_l$site <- s
    # 
    # if(j == 1){
    #   mcl_results <- mcl_site
    #   z_l_all <- z_l
    # } else{
    #   mcl_results <- bind_rows(mcl_results, mcl_site)
    #   z_l_all <- bind_rows(z_l_all, z_l)
    # }
    
    # Semivariogram analysis --------------------------------------------------
  #   print(paste("Doing semivariogram analysis for site", s))
  #   # Turn data in to a spatial dataframe for analysis
  #   xy <- xyzSub[, c(1, 2)]
  #   spdf <- SpatialPointsDataFrame(coords = xy,
  #                                  data = xyzSub,
  #                                  proj4string =
  #                                    CRS("+proj=utm +zone=15 +datum=WGS84"))
  #   
  #   # Subsample data for faster computation
  #   spSample <- spdf[sample(1:length(spdf), 10000), ]
  #   
  #   # basic variogram
  #   vario <- variogram(z~1, 
  #                      data = spSample, 
  #                      width = 0.5)
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
  #   v_fit <- fit.variogram(vario, vgm(c("Exp", 
  #                                       "Gau", 
  #                                       "Sph")))
  #   v_fit$site <- s
  #   # Sill of data
  #   zsill <- var(spSample$z)
  #   
  #   # Plot data
  #   semi <- ggplot(data = vario, aes(x = dist, y = gamma)) +
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
  #   # Collect all semivariogram data in one place
  #   semi_site <- data.frame(dist = vario$dist,
  #                           gam = vario$gamma,
  #                           site = s,
  #                           sill = zsill)
  #   if(j == 1){
  #     semi_results <- semi_site
  #     semi_fits <- v_fit
  #   } else{
  #     semi_results <- bind_rows(semi_results, semi_site)
  #     semi_fits <- bind_rows(semi_fits, v_fit)
  #   }
  }
  # print("Writing cluster analysis to .csv")
  # write.csv(mcl_results, file = paste0("mcl_results_", 
  #                                      trend, ".csv"))
  # 
  # print("Writing semivariogram analysis to .csv")
  # write.csv(semi_results, file = paste0("semi_results_", 
  #                                       trend, ".csv"))
  # write.csv(semi_fits, file = paste0("semi_fit_results_", 
  #                                    trend, ".csv"))
  # 
  # print("Plotting semivariogram analysis")
  # Remove a bunch of data
  rm(bimdl, dmcl, hu, mcl, mcl_bic, mcl_results, mcl_site,
     mcl_tidy, p_hum, semi, semi_fits, semi_site, spdf,
     spSample, trimdl, unimdl, v_fit, vario, vario_dir, xy, xyz,
     z_l, xyzSample)
  # Write the z_l_data to disc
  # write_rds(z_l_all,
  #           path = paste0("hummock_elevations_",
  #                         trend))
  
  # Add a column for plotting about site type
  # semi_results <- semi_results %>%
  #   mutate(type = str_sub(site, start = 1L, end = 1L))
  # # Plot all data on one plot for comparison
  # semi_all <- ggplot(data = semi_results, 
  #                    aes(x = dist, 
  #                        y = gam,
  #                        colour = site,
  #                        shape = type)) +
  #   scale_shape_manual(breaks = c("D", "L", "T"),
  #                      values = c(16, 1, 0),
  #                      guide = FALSE) + 
  #   scale_colour_viridis_d(name = "Site") + 
  #   geom_point(size = 1.8) + theme_bw() + 
  #   theme(legend.position = "bottom",
  #         axis.title.x = element_text(face = "bold", 
  #                                     vjust = 0.6), 
  #         axis.title.y = element_text(face = "bold", 
  #                                     vjust = 0.6),
  #         panel.grid.major.x = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         panel.grid.major.y = element_blank()) +
  #   ylab("Semivariance") +
  #   xlab("Distance (m)")
  # 
  # ggsave(plot = semi_all, 
  #        filename = paste(trend,
  #                         "all_semivariograms.tiff",
  #                         sep = "_"),
  #        device = "tiff",
  #        width = 4, height = 3, 
  #        units = "in")
  # # Remove data
  # rm(semi_all, semi_results)
  
  # Plot data colored by hummock or hollow
  # p_hum_all <- ggplot() +
  #   geom_density(data = z_l_all,
  #                aes(x = z,
  #                    y = (..scaled..) * n,
  #                    fill = as.factor(hum)),
  #                alpha = 0.3) + 
  #   scale_fill_viridis(breaks = c(0, 1),
  #                      labels = c("Hollow",
  #                                 "Hummock"),
  #                      discrete = TRUE,
  #                      option = "cividis") +
  #   theme_bw() +
  #   # geom_text(x = 0.15,
  #   #           y = 0.65,
  #   #           aes(label = paste("list(Delta*bar(z) ==", 
  #   #                             round(dif, digits=2), ")")),
  #   #           show.legend = FALSE,
  #   #           parse = TRUE) +
  #   theme(legend.position = c(0.85, 0.2),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         legend.title = element_blank()
  #   ) +
  #   facet_wrap(~site) +
  #   ylab("Count") +
  #   xlab("Relative Elevation (m)")
  # 
  # ggsave(plot = p_hum_all,
  #        filename = paste("all_sites", 
  #                         trend,
  #                         "histogram_sample.tiff",
  #                         sep = "_"),
  #        device = "tiff",
  #        width = 6, height = 4,
  #        units = "in")
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

data_path <- "RDS_files/"
files <- dir(data_path, pattern = "*")
files <- files[!(files %in% "Samples")]

# Read in all files and merge
df <- tibble(filename = files,
             site = str_sub(filename, 1, 2),
             trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% nd,
                trend == "end") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "no")

df <- tibble(filename = files,
                site = str_sub(filename, 1, 2),
                trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% l,
                trend == "ear") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "linear") %>%
  bind_rows(df)

df <- tibble(filename = files,
             site = str_sub(filename, 1, 2),
             trend = str_sub(filename, -3, -1)) %>% 
  dplyr::filter(site %in% q,
                trend == "uad") %>%
  mutate(data = map(filename,
                    ~ read_rds(file.path(data_path, .)))) %>%
  dplyr::select(-filename) %>%
  unnest() %>%
  dplyr::select(-trend) %>%
  mutate(trend = "quad") %>%
  bind_rows(df)

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
  mutate(zmaxn = ifelse(site %in% c("B1", "B3", "B6"),
                        z80n,
                        zmaxn)) %>%
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

# Subsample data for faster processing
dfSam <- df %>%
  group_by(site) %>%
  sample_n(size = 10000)
rm(df)
# Plot data binary
p_hum_best_d <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "D"),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 8500),
                     breaks = seq(0, 8000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.43,
            y = 7900,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.53,
            y = 6850,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.48,
            y = 5700,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
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
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 10000),
                     breaks = seq(0, 10000, 2500)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.43,
            y = 9400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.53,
            y = 8200,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.48,
            y = 6900,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
theme(legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      axis.text.y = element_blank()
) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_t <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "T"),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 8000),
                     breaks = seq(0, 8000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.43,
            y = 7400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.53,
            y = 6500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.48,
            y = 5450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

legend <- cowplot::get_legend(p_hum_best_d)

p_hum_best2 <- ggdraw() +
  draw_plot(p_hum_best_da + rremove("x.text") + rremove("x.title"), 
            x = 0, y = 0.7, width = 1, height = 0.3) +
  draw_plot(p_hum_best_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.4, width = 0.76, height = 0.3) +
  draw_plot(p_hum_best_t,
            x = 0, y = 0, width = 0.76, height = 0.4) + 
  draw_grob(legend,
            0.82, 0, .3/3.3, 0.5)

ggsave(plot = p_hum_best2,
       filename = "all_sites_densities_binary.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")

# Plot data with hollows classified
p_hum_best_d2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "D",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 7000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "D"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.justification = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_da2 <- p_hum_best_d2 + theme(legend.position = "none")

p_hum_best_l2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "L",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 6000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "L"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

p_hum_best_t2 <- ggplot() +
  geom_density(data = dplyr::filter(dfSam, type == "T",
                                    !is.na(hum2)),
               aes(x = z,
                   y = (..scaled..) * n,
                   fill = as.factor(hum2)),
               alpha = 0.3) + 
  scale_fill_viridis(breaks = c(0, 1),
                     labels = c("Hollow",
                                "Hummock"),
                     discrete = TRUE,
                     option = "cividis") +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = seq(0, 7000, 2000)) +
  scale_x_continuous(limits = c(-0.2, 0.65),
                     breaks = seq(-0.2, 0.6, 0.2)) +
  theme_bw() +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.43,
            y = 6400,
            aes(label = paste("list(Delta*bar(z) ==",
                              round(zavg, digits = 2),
                              "%+-%", 
                              round(zsd, digits = 2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.53,
            y = 5500,
            aes(label = paste("list(n ==",
                              no, ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  geom_text(data = dplyr::filter(df_sum, type == "T"),
            x = 0.48,
            y = 4450,
            aes(label = paste("list(A[ratio] ==",
                              round(aratio, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE,
            size = 2.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) +
  facet_wrap(~site, ncol = 4, scales = "free_y") +
  ylab("Density") +
  xlab("Relative Elevation (m)")

legend2 <- cowplot::get_legend(p_hum_best_d)

p_hum_best3 <- ggdraw() +
  draw_plot(p_hum_best_da2 + rremove("x.text") + rremove("x.title"), 
            x = 0, y = 0.7, width = 1, height = 0.3) +
  draw_plot(p_hum_best_l2 + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.4, width = 0.76, height = 0.3) +
  draw_plot(p_hum_best_t2,
            x = 0, y = 0, width = 0.76, height = 0.4) + 
  draw_grob(legend2,
            0.82, 0, .3/3.3, 0.5)

ggsave(plot = p_hum_best3,
       filename = "all_sites_densities_hollow_classified.tiff",
       device = "tiff",
       width = 6, height = 4,
       units = "in")
