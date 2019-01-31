# 
# Author: Jake Diamond
# Purpose: Analysis of hummocks, point processes
# Date: December 5, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(MASS)
library(spatstat)
library(tidyverse)
library(broom)
library(scales)
library(sf)


# Get hummock data (x, y, area, vol, perim)
# df <- read.csv("Lidar/delineate_all_sites.csv")
df <- read.csv("hummock_stats.csv")

# # Load hummock raster
d2 <- raster("C:/Users/diamo/Dropbox/Projects/EAB/delineation_code/output/D2_hum.tif")
e <- extent(d2)
d2 <- st_as_sf(d2)


# Clean data
df <- df %>%
  dplyr::select(-X) %>%
  mutate(site = str_sub(id, 1, 2))

# Total hummock areas without cleaning
areas <- df %>%
  group_by(site) %>%
  summarize(a_hums = sum(area_poly),
            a_site = mean(site_area),
            ratio = a_hums / a_site)

# Hummock spatial analysis ------------------------------------------------
df_n <- df %>%
  group_by(site) %>%
  nest()

df_d2 <- dplyr::filter(df, 
                       site == "B1",
                       area_poly > 0.3,
                       PA_poly < 2)
ranges <- owin(xrange = c(quantile(df_d2$centroid_x, 0.05),
                          quantile(df_d2$centroid_x, 0.95)),
               yrange = c(quantile(df_d2$centroid_y, 0.05),
                          quantile(df_d2$centroid_y, 0.95))
               )
bounds <- boundingbox(ranges)

d2_ppp <- ppp(x = df_d2$centroid_x,
              y = df_d2$centroid_y,
              window = bounds,
              marks = df_d2$area_poly)
unitname(d2_ppp) <- c("meter", "meter")
# d2_ppp <- rotate.ppp(d2_ppp, angle = pi/5)

# Plot marked spatial point process
plot(d2_ppp)

# Calculate nearest neighbor distances
nn <- nndist(d2_ppp)
nn_avg <- mean(nn)
hist(nn, breaks = 20)

# Expected nearest neighbor distance between randomly dispersed
# pairs (0.5 / sqrt(n/D)), D is domain size and n is no. of hums
n_pts <- npoints(d2_ppp)
D <- area.owin(bounds)
nn_exp <- 0.5 / sqrt(n_pts / D)

# Standard error of the nearest neighbor distribution
se_nn <- 0.26 / sqrt(n_pts^2 / D)

# Ratios of observed to expected nearest neighbor
# Values greater than 1 = overdispersion
# Values below 1 = clustering
nn_ratio <- nn / nn_exp
hist(nn_ratio, breaks = 20)
nn_ratio_avg <- nn_avg / nn_exp

# z-scores to evaluate significance of overdispersion/clustering
nn_z <- (nn - nn_exp) / se_nn
hist(nn_z)
nn_z_avg <- (nn_avg - nn_exp) / se_nn
pval_avg <- 2 * pnorm(-abs(nn_z_avg))
p_val <- 2 * pnorm(-abs(nn_z))
hist(p_val)

marks(d2_ppp) <- nndist(d2_ppp)
plot(d2_ppp, markscale=0.5)
plot(Gest(d2_ppp))

# KS test against complete spatial randomness
KS <- cdf.test(unmark(d2_ppp), "y")
KS

# Test against Poisson homogeneity
lambda <- function(x, y){
  100 * (x + y)
}
plot(rpoispp(lambda))

# Homogenous poisson model fit to data
ppm(unmark(d2_ppp), ~1)

# Inhomogenous poisson model with intensity that is log-linear in the cartestian coordinates
fit <- ppm(unmark(d2_ppp), ~x + y)
fit
plot(fit, how = "image")
M.fit <- quadrat.test(fit,
                      nx = 5,
                      ny = 5)
M.fit
plot(unmark(d2_ppp), pch = ".")
plot(M.fit, add = T,  col = "red")
# Analysis of deviance for poisson
fitnull <- update(fit, ~1)
anova(fitnull, fit, test = "Chi")
AIC(fit)
AIC(fitnull)

X <- rmh(fit)
plot(X)
plot(unmark(d2_ppp))


# distance between points
pairdist(d2_ppp)
nndist(d2_ppp)
emp <- distmap(d2_ppp)
plot(emp, main = "Empty space distances")
plot(d2_ppp, add = TRUE)


# 
Fc <- Fest(d2_ppp)
Fc


# Nearest neighbor comparison
Gc <- Gest(d2_ppp)
Gc
par(pty = "s")
plot(Gc)

# Pairwise distances and K function
Gc.k <- Kest(d2_ppp)
Gc.k
par(pty = "s")
plot(Gc.k)

# 
L <- Lest(d2_ppp)
plot(L, main = "L function")

# Pair correlation = probability of observing a pair of points separated by distance r, divided by corresponding probability for a Poisson process, g(r) = 1 is complete randomness, g(r) > 1 is clustering or attraction, g(r) < 1 is inhibition or regularity
plot(pcf(d2_ppp))

# J function, Values J(r) > 1 suggest regularity, and J(r) < 1 suggest clustering
J <- Jest(d2_ppp)
plot(J)

# Plot all functions
plot(allstats(unmark(d2_ppp)))

# Envelopes
E <- envelope(d2_ppp, Lest, nsim = 19, rank = 1, global = TRUE)
E
plot(E)


# Inhomogeneous process
lam <- predict(fit, locations = d2_ppp)
Ki <- Kinhom(d2_ppp, lam)
plot(Ki)







# Hummock distribution analysis -------------------------------------------
# Get data in descending rank order
df_h <- df %>%
  dplyr::select(-(1:10),-(12:14),-(16:17), -21) %>%
  gather(key = "measure.type", value = "measure",
         -site) %>%
  group_by(site, measure.type) %>%
  mutate(rank = cume_dist(desc(measure)))

# Model for exponential estimate model
e_mod <- function(data) {
  lm(rank ~ exp(measure), data = data)
}

# Fit exponential line to area, perim, vol
mod <- df_h %>%
  group_by(site, measure.type) %>%
  nest() %>%
  mutate(model = map(data, e_mod),
         tidied = map(model, tidy)) %>% 
  unnest(tidied, .drop = TRUE) %>%
  ungroup()

df_h$measure.type <- factor(df_h$measure.type,
                            levels = c("perimeter",
                                       "area",
                                       "volume"))
levels(df_h$measure.type) <- c("Perimeter~(m)",
                               "Area~(m^{2})",
                               "Volume~(m^{3})")

# Plot
p_rank <- ggplot(data = df_h,
                aes(x = measure,
                    y = rank,
                    color = site)) +
  geom_point(size = 3) +
  # scale_shape_manual(name = "Site",
  #                    values = c(1, 16)) +
  # stat_smooth(formula = y ~ exp(x)) +
  facet_wrap(~measure.type,
             labeller = label_parsed) +
  theme_bw() +
  scale_x_log10(
    labels = trans_format('log10', math_format(10 ^ .x))
    ,
    limits = c(10 ^ -1.2, 10 ^ 1.2),
    breaks = c(10 ^ -1, 10 ^ -0, 10 ^ 1)
    ) +
  scale_y_log10(
    labels = trans_format('log10', math_format(10 ^ .x))
    ,
    limits = c(10 ^ -2.2, 10 ^ 0.2),
    breaks = c(10 ^ -2, 10 ^ -1, 10 ^ 0)
    ) +
  annotation_logticks(base = 10, sides = "tlbr") + 
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 18, colour = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = c(0.12, 0.17),
    legend.text = element_text(size = 16),
    legend.background = element_rect(
      fill = "gray90",
      size = 1,
      linetype = "solid",
      colour = "black")) + 
  xlab("") +
  ylab(expression("P (" * X>=x * ")"))

p_rank

ggsave(p_rank, filename = "Hummock_exponential.tiff",
       device = "tiff",
       dpi = 600,
       width = 10,
       height = 6)



r <- raster("C:/Users/diamo/Dropbox/Projects/EAB/delineation_code/output/D1_hum_ho_binary.tif")
plot(r)
df <- rasterToPoints(r)
df <- as.data.frame(df)
colnames(df)[3] <- "hu_ho"
mod_binom = glm(hu_ho ~ x+y, data = df, family = binomial)  
summary(mod_binom)

library(spectral)
x <- sort(unique(df$x))
y <- sort(unique(df$y))

n = length(unique(df$x))
m = length(unique(df$y))
mat = matrix(0, nrow=m, ncol=n)
for (i in nrow(df)){
  mat[df[i,"x"], df[i,"y"]] = df[i,"hu_ho"]
}
library(kzfs)
# Raw periodogram (2d spectrum)
rp <- kzp2(mat)
rp_sum <- kzp2
fy <- rp$freq.y
fx <- rp$freq.x
rp <- rp$kzp2d
# smoothing 2D spectrum 2 times
sp <- smooth.kzp2(rp,0.01,k=2)
par(mfrow=c(2,1), cex=0.5)
persp(x=fx, y=fy, z=rp, expand =0.5,
      main = "Raw 2D KZ Periodogram", ltheta=40, shade=0.75,
      theta=-30, phi=15, zlab="",xlab="x", ylab="y",
      ticktype="detailed", col="red")
persp(x=fx, y=fy, z=sp, expand =0.5,
      main = "Smoothed 2D KZ Periodogram", ltheta=40, shade=0.75,
      theta=-30, phi=25, zlab="",xlab="x", ylab="y",
      ticktype="detailed", col="red")
par(mfrow=c(1,1), cex=1)
kzp2.summary(sp) # direction & frequency




FT <- spec.fft(x = round(abs(min(x)) + x, 2), 
               y = round(abs(min(y)) + y, 2),
               z = mat)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x,
             y = y,
             z = m,
             main = "Propagating Wave")
plot(FT)

# x <- round((sort(unique(df$x)) + abs(min(df$x))) * 1000, 0)
# y <- round((sort(unique(df$y)) + abs(min(df$y))) * 1000, 0)
# m <- matrix(NA, length(x) + 1, length(y) + 1)
# m[cbind(x, y)] <- df$hu_ho
# 
# mat <- matrix(NA, nrow=max(df[[1]])+1, ncol=max(df[[2]])+1 ) 
# mat[data.matrix(df[,1:2])] <- df[,3] 
# mat[cbind(x, y)] <- df$hu_ho
# mat2 <- data.matrix(df)
# 
# mtx <- matrix(NA, nrow=length(unique(df$x)) + 2, ncol=length(unique(df$y)) + 2 )
# mtx[cbind(order(df$x), order(df$y))] <- df$hu_ho




# Test with data from package

x <- seq(0, 1, length.out = 1e2)
y <- seq(0, 1, length.out = 1e2)
# calculate the data
m <- matrix(0, length(x), length(y))
for (i in 1:length(x))
  for (j in 1:length(y))
    m[i, j] <- sin(4 * 2 * pi * x[i] + 10 * 2 * pi * y[j])
# calculate the spectrum
FT <- spec.fft(x = x, y = y, z = m)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x,
             y = y,
             z = m,
             main = "Propagating Wave")
plot(
  FT,
  main = "2D Spectrum",
  palette = "wb"
  ,
  xlim = c(-20, 20),
  ylim = c(-20, 20),
  zlim = c(0, 0.51)
  ,
  xlab = "fx",
  ylab = "fy",
  zlab = "A",
  ndz = 3,
  z.adj = c(0, 0.5)
  ,
  z.cex = 1
)


# Test with data
x <- sort(unique(df$x))
y <- sort(unique(df$y))

n = length(unique(df$x))
m = length(unique(df$y))
mat = matrix(0, nrow=m, ncol=n)
for (i in nrow(df)){
  mat[df[i,"x"], df[i,"y"]] = df[i,"hu_ho"]
}


x <- sort(unique(round(abs(min(df$x)) + df$x, 2)))
y <- sort(unique(round(abs(min(df$y)) + df$y, 2)))

x2 <- seq(0, 51.3, length.out = 360)
y2 <- seq(0, 54, length.out = 180)

mat2 <- as.matrix(r)
# calculate the data
m2 <- mat
# smooth data for irregular data
library(fields)
mat_smooth <- image.smooth(mat2)
mat3 <- mat_smooth$z
mat3 <- mat3[complete.cases(mat3),complete.cases(t(mat3))]

# calculate the spectrum
FT <- spec.fft(x = x2, y = y2, z = mat3)
# plot
par(mfrow = c(2, 1))
rasterImage2(x = x2,
             y = y2,
             z = mat2,
             main = "Propagating Wave")
plot(
  FT,
  main = "2D Spectrum",
  palette = "wb"
  ,
  # xlim = c(-20, 20),
  # ylim = c(-20, 20),
  # zlim = c(0, 0.51)
  # ,
  xlab = "fx",
  ylab = "fy",
  zlab = "A",
  # ndz = 3,
  # z.adj = c(0, 0.5)
  # ,
  z.cex = 1
)
dev.off()
