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
library(tidyverse)
library(mixtools)
library(MASS)
library(gstat)

# Get all filenames
filenames <- paste("Lidar/detrended/", 
                   list.files("Lidar/detrended"), 
                   sep = "")
# Bimodal analysis --------------------------------------------------------
for (i in 1:length(filenames)){
  rm(p_b, p_u)
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
  # Convert raster to points
  xyz <- rasterToPoints(r)
  xyz <- as.data.frame(xyz)
  colnames(xyz)[3] <- "z"
  xyz <- dplyr::filter(xyz, z < quantile(z, 0.9))
  xy <- xyz[, c(1, 2)]
  spdf <- SpatialPointsDataFrame(coords = xy, 
                                 data = xyz,
                                 proj4string = 
                                   CRS("+proj=utm +zone=15 +datum=WGS84"))
  spSample <- spdf[sample(1:length(spdf), 100000), ]
  
  # Remove data for memory
  rm(r)
  
  # Get a bimodal model of two normal distributions
  mixmdl <- normalmixEM(na.omit(spSample$z),
                        k = 2
                        )
  # Get a fit for a unimodal normal distribution
  unimdl <- fitdistr(spSample$z, "normal")
  para <- unimdl$estimate
  # Get log-likelihoods for both uni and bimodal distributions
  uni <- unimdl$loglik
  bi <- mixmdl$loglik
  # Calculate the deviance between them
  D <- bi - uni
  # Determinen  if significantly different
  p <- 1 - pchisq(D, df = 3)
  # Remove some data for memory
  rm(spdf, spSample)
  
  # Plot bimodal if p <0.01, else unimodal
  if(p < 0.01){
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
      stat_function(fun = plot_mix_comps, aes(colour = "1"),
                    args = list(mixmdl$mu[1], mixmdl$sigma[1],
                                lam = mixmdl$lambda[1]), lwd = 1.5) +
      stat_function(fun = plot_mix_comps, aes(colour = "2"),
                    args = list(mixmdl$mu[2], mixmdl$sigma[2],
                                lam = mixmdl$lambda[2]), lwd = 1.5) +
      geom_segment(aes(x = mixmdl$mu[1],
                   xend = mixmdl$mu[1],
                   y = 0,
                   yend = max(density(xyz$z)$y)),
                   colour = "blue",
                   linetype = "dashed",
                   size = 1.2
                   ) +
      geom_segment(aes(x = mixmdl$mu[2],
                   xend = mixmdl$mu[2],
                   y = 0,
                   yend = max(density(xyz$z)$y)),
                   colour = "green",
                   linetype = "dashed",
                   size = 1.2
                   ) +
      geom_text(aes(x = mixmdl$mu[1] + 0.08,
                    y = max(density(xyz$z)$y) + 0.05,
                    label = round(mixmdl$mu[1], 2)),
                color = "blue"
                ) + 
      geom_text(aes(x = mixmdl$mu[2] + 0.08,
                    y = max(density(xyz$z)$y) + 0.05,
                    label = round(mixmdl$mu[2], 2)),
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
                                        vjust = 0.6,
                                        size = 28), 
            axis.title.y = element_text(face = "bold", 
                                        vjust = 0.6,
                                        size = 28),
            axis.text.x = element_text(size = 24),
            axis.text.y = element_text(size = 24),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            ) +
      ylab("Density") +
      xlab("Relative Elevation (m)")
    
    p_b
    ggsave(plot = p_b, filename = paste0(s, "_Bimodal_detrended.tiff"),
           device = "tiff",
           width = 8, height = 6, 
           units = "in")
  } else {
    # Function to plot normal distributions
    plot_comps <- function(x, mu, sigma, lam) {
      dnorm(x, mu, sigma)
    }
    # Plotting unimodal
    p_u <- ggplot(data = xyz, aes(x = z,
                                fill = ..x..)) +
      geom_histogram(aes(z, 
                         ..density..), 
                     bins = 100) +
      stat_function(fun = plot_comps,
                    args = list(unimdl$estimate[1], unimdl$estimate[2]), 
                    lwd = 1.5,
                    color = "black") +
      scale_fill_gradientn(colours = c("blue", 
                                       "green", 
                                       "yellow",
                                       "red")) +
      theme_bw() + 
      theme(legend.position = "none",
            axis.title.x = element_text(face = "bold", 
                                        vjust = 0.6,
                                        size = 28), 
            axis.title.y = element_text(face = "bold", 
                                        vjust = 0.6,
                                        size = 28),
            axis.text.x = element_text(size = 24),
            axis.text.y = element_text(size = 24),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            ) +
      ylab("Density") +
      xlab("Relative Elevation (m)")
    
    p_u
    ggsave(plot = p_u, filename = paste0(s, "_Unimodal_detrended.tiff"),
           device = "tiff",
           width = 8, height = 6, 
           units = "in")
  }
}


# Semivariogram analysis --------------------------------------------------
vario <- variogram(D1_1cm ~ 1,
                   data = spSample)
zsill <- var(spdf$D1_1cm)
plot(vario$dist, vario$gamma, 
     # xlim=c(0,8000),
     # ylim = c(0, zsill + 1),
     xlab="Distance (m)",
     ylab="Semivariogram")
abline(h = zsill)
text(2, 
     zsill + 0.001, 
     paste("Sill =", round(zsill, 3)))





locminlocations <- which(diff(sign(diff(zvario$gamma)))==-2)+1
locmin <- zvario$dist[locminlocations]
locmin <- locmin[-c(1:6, 8:11, 13:19, 21:22, 24:25)]


zvm <- fit.variogram(gstat_variogram, model = vgm(c("Exp", "Sph")),
                     fit.sills = FALSE)
curve(0.02117 * (1 - exp(-x / 2.638637)), add = TRUE)

semi <- ggplot(data = zvario, aes(x = dist, y = gamma)) +
  geom_point(size = 1.8) + theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(face = "bold", 
                                    vjust = 0.6,
                                    size = 28), 
        axis.title.y = element_text(face = "bold", 
                                    vjust = 0.6,
                                    size = 28),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  geom_hline(yintercept =  zsill, 
             size = 1.4,
             linetype = "dashed") + 
  geom_vline(xintercept = locmin,
             size = 1,
             linetype = "dashed",
             color = "dark grey") +
  annotate("text", x = 2, y = 0.025, 
           label = paste("Sill =", 
                         round(zsill, 3)),
           size = 8) +
  scale_x_continuous(breaks = seq(0, 15, 2)) +
  scale_y_continuous(breaks = seq(0, 0.03, 0.005)) +
  
  ylab("Semivariogram") +
  xlab("Distance (m)")
ggsave(plot = semi, filename = "D2_semivariogram.png",
       width = 8, height = 6, 
       units = "in")