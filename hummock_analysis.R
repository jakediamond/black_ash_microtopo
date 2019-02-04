# 
# Author: Jake Diamond
# Purpose: Analysis of hummocks, point processes
# Date: December 5, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(data.table)
library(raster)
library(MASS)
library(spatstat)
library(tidyverse)
library(broom)
library(scales)
library(sf)


# Get hummock data (x, y, area, vol, perim)
# df <- read.csv("Lidar/delineate_all_sites.csv")
df <- read.csv("Lidar/hummock_stats_ext6.csv") %>%
  dplyr::select(id, site_area) %>%
  mutate(site = str_sub(id, 1, 2)) %>%
  select(site, site_area) %>%
  distinct() %>%
  right_join(read.csv("hummock_stats_clean.csv"))

# Total hummock areas without cleaning
areas <- df %>%
  group_by(site) %>%
  summarize(a_hums = sum(area),
            a_site = mean(site_area),
            ratio = a_hums / a_site)

# Minimum area for hummocks to consider
min_hum <- 0.1

# Total hummock areas with cleaning
areas_clean <- df %>%
  group_by(site) %>%
  dplyr::filter(area > min_hum) %>%
  summarize(a_hums = sum(area),
            a_site = mean(site_area),
            ratio = a_hums / a_site)



# df_hums <- fread("Lidar/D1_hum_unnormalized_ext6_cleaned.txt")
# colnames(df_hums) <- c("x", "y", "z", "id", "class")
# 
# # Plot all data
# hums_p <- ggplot(data = df_hums,
#                  aes(x = x,
#                      y = y)) + 
#   # facet_wrap(~site) + 
#   geom_polygon(aes(fill = id, group = id)) + 
#   labs(x = "Long.", y = "Lat.", 
#        title = "Samples from the predictive process")
# hums_p
# Hummock spatial analysis ------------------------------------------------
# Nest data by site to apply the spatial analysis function to
df_n <- df %>%
  group_by(site) %>%
  nest()

# Function to do hummock analysis
# Includes plot of marked process, nearest neighbors results
# and plot. Also plots pair-correlation function and G_est
hum_spat_fun <- function(data){
  site <- unique(data$site)
  # Get polygon window of data (needs to be in 
  # anitclockwise order)
  poly <- rev(chull(x = data$xmean, y = data$ymean))
  polyrange <- owin(poly = list(x = data$xmean[poly],
                                y = data$ymean[poly]))
  # Create the marked point process for analysis
  pp <- ppp(x = data$xmean,
            y = data$ymean,
            window = polyrange,
            marks = data$area)
  unitname(pp) <- c("meter", "meter")
  # Calculate nearest neighbor distances
  nn <- nndist(pp)
  nn_avg <- mean(nn)
  # hist(nn, breaks = 20)
  
  # Expected nearest neighbor distance between randomly dispersed
  # pairs (0.5 / sqrt(n/D)), D is domain size and n is no. of hums
  n_pts <- npoints(pp)
  D <- area.owin(polyrange)
  nn_exp <- 0.5 / sqrt(n_pts / D)
  
  # Standard error of the nearest neighbor distribution
  se_nn <- 0.26 / sqrt(n_pts^2 / D)
  
  # Ratios of observed to expected nearest neighbor
  # Values greater than 1 = overdispersion
  # Values below 1 = clustering
  nn_ratio <- nn / nn_exp
  nn_ratio_avg <- nn_avg / nn_exp
  
  # z-scores to evaluate significance of overdispersion/clustering
  nn_z <- (nn - nn_exp) / se_nn
  nn_z_avg <- (nn_avg - nn_exp) / se_nn
  pval_avg <- 2 * pnorm(-abs(nn_z_avg))
  p_val <- 2 * pnorm(-abs(nn_z))
  
  
  marks(pp) <- nndist(pp)
  plot(pp, markscale=0.5)
  plot(Gest(pp))
  
  # distance between points
  emp <- distmap(pp)
  plot(emp, main = "Empty space distances")
  plot(pp, add = TRUE)
  
  # Pair correlation = probability of observing a pair of points separated by distance r, divided by corresponding probability for a Poisson process, g(r) = 1 is complete randomness, g(r) > 1 is clustering or attraction, g(r) < 1 is inhibition or regularity
  plot(pcf(pp))
  data.frame(nn = nn,
             nn_avg = nn_avg,
             n_pts = n_pts,
             D = D,
             nn_exp = nn_exp,
             se_nn = se_nn,
             nn_ratio = nn_ratio,
             nn_ratio_avg = nn_ratio_avg,
             pval_avg = pval_avg,
             pval = p_val
  )
}

test <- df_n %>%
  mutate(pps = map(data, hum_spat_fun))
test2 <- unnest(test, pps, .drop = TRUE)



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
E <- envelope(d2_ppp, Gest, nsim = 19, rank = 1, global = TRUE)
E
plot(E)

# Hummock distribution analysis -------------------------------------------
# Get data in descending rank order
df_h <- df %>%
  dplyr::select(site, area, area_poly,
                vol, perim_poly) %>%
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
                            levels = c("perim_poly",
                                       "area",
                                       "area_poly",
                                       "vol"))
levels(df_h$measure.type) <- c("Perimeter~(m)",
                               "Area~(m^{2})",
                               "AreaP~(m^{2})",
                               "Volume~(m^{3})")
# don't plot these sites
bad <- c("B1", "B3", "B6")

# Plot
p_rank <- ggplot(data = filter(df_h,
                               !(site %in% bad)),
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




# KS test against complete spatial randomness
KS <- cdf.test(unmark(d2_ppp), "x")
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

# Inhomogeneous process
lam <- predict(fit, locations = d2_ppp)
Ki <- Kinhom(d2_ppp, lam)
plot(Ki)