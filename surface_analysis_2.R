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
library(ecogen)

# Get raster files from point cloud 
dfr <- raster("Lidar/Rasters/D1_1cm.tif")
projection(dfr)
plot(dfr)
xyz <- rasterToPoints(dfr)
xyz <- as.data.frame(xyz)










# xyz$D1_1cm <- ifelse(xyz$D1_1cm == min(xyz$D1_1cm),
#                      NA,
#                      xyz$D1_1cm)
# xyz <- na.omit(xyz)
xy <- xyz[, c(1, 2)]
spdf <- SpatialPointsDataFrame(coords = xy, 
                               data = xyz,
                               proj4string = 
                                 CRS("+proj=utm +zone=15 +datum=WGS84"))



plot(hist(xyz$D1_1cm))

point_data <- as(dfr, 'SpatialPointsDataFrame')
point_data$D1_1cm <- ifelse(point_data$D1_1cm == min(point_data$D1_1cm),
                            NA,
                            point_data$D1_1cm)
point_data <- (point_data$)
gstat_variogram <- variogram(D1_1cm ~ 1, ,
                             data = spdf,
                             width = 2)





df <- read.table("Lidar/Clouds/D1_1cm.asc")
df <- select(df, 1:3) %>%
  rename(x = V1, y = V2, z = V3)
df_z.d <- select(df, 1:2, 7) %>%
  rename(x = V1, y = V2, z.d = V10)
e <- extent(df_z[,1:2])
plot(df)
# Detrending
# Fit a linear model with the x and y coordinates as predictors
fit.test <- lm(z ~ x + y, data=df_z)

# Extract the residuals
lm.residuals <- fit.test$residuals  

# To fit polynomial of degree two or higher, you must pass a matrix with the xy coordinates through poly()
# isolate the columns with the x and y coordinates and covert to matrix
X.Y <- as.matrix(df_z[,c(1,2)]) 
# Create a vector with the elevation data
z <- df_z[, 3]  

#Recreate the data frame with the coordinates and the residuals (instead of the raw elevation)

x.y <- df_z[,1:2]  #isolate the coordinate columns

x.y$residuals <- lm.residuals #create column of residuals
x.y$quad  <- quad.residuals
x.y.r <- x.y   #rename the data frame

head(x.y.r)  #check to make sure it worked
plot(hist(df_z$z, breaks= 1000), freq = FALSE)
plot(hist(df$V10, breaks= 1000), freq = FALSE)
plot(hist(x.y.r$residuals, breaks= 1000), freq = FALSE)
plot(hist(x.y.r$quad, breaks= 1000), freq = FALSE)

#convert back to RasterLayer and restore the projection and datum
df.detrend <- rasterFromXYZ(x.y.r, 
                            res = c(0.01, 0.01)
                            )

#convert back to SpatialGridDataFrame

alex.15.detrended <- as(alex.15.detrend, ‘SpatialGridDataFrame’)

#convert to im object in spatstat

alex.15.d <- as.im(alex.15.detrended)




coordinates(df_z) <- ~x+y
r <- raster(ext = extent(df_z), resolution = 1)
r_z <- rasterize(df_z, r, df_z$z, fun = mean)
# r_z.d <- rasterize(df_z.d[,1:2], r, df_z.d[,3], fun = mean)
plot(r_z)
# plot(r_z.d)
coordinates(x.y.r) <- ~x+y
r <- raster(ext = extent(x.y.r), resolution = 1)
r_z <- rasterize(x.y.r, r, x.y.r$residuals, fun = mean)
# r_z.d <- rasterize(df_z.d[,1:2], r, df_z.d[,3], fun = mean)
plot(r_z)



points <- SpatialPoints(y[,1:2], y[,3])
# Read in xy coordinates
df_xy <- data.frame(x = c(1, 2, 3), 
                    y = c(1, 2, 3)
)
# Extract raster values from xy coordinates
xy_r_values <- extract(r3, 
                       SpatialPoints(df_xy), 
                       sp = TRUE
)@data 

# 2) you want to analyze the distributions in what way? If you just want the raw numbers to throw in a model of the size distribution just grab the Z values (raw) or the scalar field (detrended) of .asc files. You would read it like this, I believe:

dist <- data.frame(Z = read.table("my_model.asc", sep = ",")[,3], Zd =  read.table("my_model.asc", sep = ",")[,7])
#from here you can run any stats you'd like to analyze the distribution


df <- read.csv("D2_10cm_vertices1m_h_v3.csv")


# Bimodal analysis
plot(hist(df_z.d$z.d, breaks= 50), freq = FALSE)
mixmdl <- normalmixEM(na.omit(xyz$D1_1cm),
                      epsilon = 100,
                      maxit = 100)
plot(mixmdl, which = 2, breaks = 72)
,
     xlab = "Elevation (m)",
     ylab = "Density")
lines(density(df_z$z, lty = 2, lwd = 2))

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

p <- ggplot(data = df, aes(x = Coord._Z,
                           fill = ..x..)) +
  geom_histogram(aes(Coord._Z, 
                     ..density..), 
                 bins = 72) +
  stat_function(fun = plot_mix_comps, aes(colour = "1"),
                args = list(mixmdl$mu[1], mixmdl$sigma[1], 
                            lam = mixmdl$lambda[1]), lwd = 1.5) +
  stat_function(fun = plot_mix_comps, aes(colour = "2"),
                args = list(mixmdl$mu[2], mixmdl$sigma[2], 
                            lam = mixmdl$lambda[2]), lwd = 1.5) +
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
        panel.grid.major.x = element_line(size = 0.8,
                                          linetype = "dashed",
                                          color = "dark grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(size = 0.8,
                                          linetype = "dashed",
                                          color = "dark grey")) +
  ylab("Density") +
  xlab("Relative Elevation (m)")
ggsave(plot = p, filename = "D2_Bimodal.png",
       width = 8, height = 6, 
       units = "in")


fit <- fitdistr(df$z, "normal")
para <- fit$estimate
uni <- fit$loglik
bi <- mixmdl$loglik
D <- 2 * (bi - uni)

pchisq(D, df=2, lower.tail=FALSE)

# Semivariogram analysis
points <- SpatialPoints(xyz[,1:2], xyz[,3])
xyspatial <- SpatialPoints(cbind(df_z$x, df_z$y))
zspatial <- data.frame(df_z$z)
spatialdata <- SpatialPointsDataFrame(xyspatial, zspatial)
spatialdata <- as(dfr, 'SpatialPointsDataFrame')
xyz <- as.data.frame(spatialdata)
xyz$D1_1cm <- ifelse(xyz$D1_1cm == min(xyz$D1_1cm),
                     NA,
                     xyz$D1_1cm)
xyz <- na.omit(xyz)
zvario <- variogram(spatialdata$D1_1cm ~ 1, 
                    locations = spatialdata,
                    width = 0.1)

zsill <- var(df$Coord._Z)
plot(zvario$dist, zvario$gamma, 
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


zvm <- fit.variogram(zvario, model = vgm("Pow", 0.02117),
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