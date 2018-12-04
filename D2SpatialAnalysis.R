#
# Purpose: To analyze D2 spatial data
# Coder: Jake Diamond
# Date: November 30, 2017
#

# Set Working directory
setwd("C:/Users/Jake/Dropbox/Projects/EAB/Data")

# Load libraries
library(dplyr)
library(mixtools)
library(MASS)
library(spatstat)
library(ggplot2)
library(raster)
library(spectral)
library(kzfs)
library(gstat)
library(ncf)
ncf.cor <- correlog(df$X..X, df$Y, df$Coord._Z,
                    increment=0.2, resamp=5)
ncf.cor
plot(x = ncf.cor$mean.of.class, y = ncf.cor$correlation)

library(SpatialPack)
data <- summary(modified.ttest(z.value,z.value,coords = xy,nclass = 21))
plot(x=data$coef[,1],y = data$coef[,4],type = "l")
# Read in elevation data
df <- read.csv("D2_10cm_vertices1m_h_v3.csv")
# df <- read.csv("D2_10cm_vertices1m_h_v3.csv")
# df <- df[df$Coord..Z < -1.3, ]

# Bimodal analysis
mixmdl <- normalmixEM(df$Coord._Z)
plot(mixmdl, which = 2, breaks = 72,
     xlab = "Elevation (m)",
     ylab = "Density")
lines(density(df$Coord._Z, lty = 2, lwd = 2))

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


fit <- fitdistr(df$Coord._Z, "normal")
para <- fit$estimate
uni <- fit$loglik
bi <- mixmdl$loglik
D <- 2 * (bi - uni)

pchisq(D, df=2, lower.tail=FALSE)



# gstat stuff
xyspatial <- SpatialPoints(cbind(df$X..X, df$Y))
zspatial <- data.frame(df$Coord._Z)
spatialdata <- SpatialPointsDataFrame(xyspatial, zspatial)
zvario <- variogram(df.Coord._Z ~ 1, 
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





# Spatstat stuff
P <- ppp(df$X..X, df$Y, c(-21.25, 23.16), c(-32.9, 20), marks = df$Coord..Z)
plot(P)
kstest(P, function(x, y, m) {
  m
})
kpp <- Kest.fft(P, 0.01)
plot(kpp)
z <- rasterFromXYZ(df, digits = 3)
z2 <- as.matrix(z$Coord..Z)
spe <- spec.fft(z = z2)
spe
df$X..X[2] - df$X..X[1]
z2 <- z2[rowSums(is.na(z2)) != ncol(z2), ]
rp <- kzp2(z2)
dx <- round(max(df$X..X) - min(df$X..X), 1) * 10 #cm
dy <- round(max(df$Y) - min(df$Y), 1) * 10 #cm
b <- expand.grid(x = 1:dx, y = 1:dy)
b$z <- df$Coord..Z
a <- array(0, c(dx, dy))
a[as.matrix(b[, 1:2])] <- b$z

dx <- 100				# x range
dy <- 120				# y range
b <- expand.grid(x=1:dx, y=1:dy)
q1 <- pi/6; f1 <- 0.2;
b$v1 <- sin(f1*2*pi*(b$x*cos(q1)+b$y*sin(q1))+100*runif(1))
q2 <- pi/4; f2 <- 0.08;
b$v2 <- sin(f2*2*pi*(b$x*cos(q2)+b$y*sin(q2))+100*runif(1))
a <- array(0,c(dx,dy))
a[as.matrix(b[,1:2])] <- b$v1 + 1.5*b$v2
a <- a + 10*matrix(rnorm(dx*dy,0,1),ncol=dy)

rp <- kzp2(a)			# raw 2D spectrum

fy <- rp$freq.y; fx <- rp$freq.x; rp <- rp$kzp2d

# smoothing 2D spectrum 2 times
sp <- smooth.kzp2(rp,0.01,k=2)	

par(mfrow=c(2,1), cex=0.5)
heatmap(x= rp, Rowv = NA, Colw = NA)


persp(x=fx, y=fy, z=rp, expand =0.5,
      main = "Smoothed 2D KZ Periodogram",
      zlab="",xlab="x", ylab="y",
      ticktype="detailed")
par(mfrow=c(1,1), cex=1)

kzp2.summary(rp)		# direction & frequency
