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
library(poweRlaw)

# Get hummock data (x, y, area, vol, perim)
# df <- read.csv("Lidar/delineate_all_sites.csv")
df <- read.csv("Lidar/hummock_stats_ext6.csv") %>%
  dplyr::select(id, site_area) %>%
  mutate(site = str_sub(id, 1, 2)) %>%
  dplyr::select(site, site_area) %>%
  distinct() %>%
  right_join(read.csv("hummock_stats_clean.csv"))
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"
df$type <- str_sub(df$site, 1, 1)
df$num <- str_sub(df$site, 2, 2)

# Minimum area for hummocks to consider
min_hum <- 0.1

# Estimate how much additional surface area is provided by 
# hummocks.
df_sa <- df %>%
  transmute(site = site,
            site_area = site_area,
            sa = z80n * perim_poly) %>%
  group_by(site) %>%
  summarize(site_area = mean(site_area),
            sa = sum(sa),
            sa_add = sa / site_area)

# Estimate how much of estimated wetland volume is taken 
# up by hummocks. Wetland volume is estimated as the mean 80th 
# quantile of hummock height in a site * wetland area
df_vol <- df %>%
  group_by(site) %>%
  summarize(site_ht = mean(z80n),
            site_vol = site_ht * mean(site_area),
            hum_vol = sum(vol),
            hum_ratio = hum_vol / site_vol)

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
  dplyr::filter(area > min_hum) %>%
  group_by(site) %>%
  nest()

# Function to do hummock analysis
# Includes plot of marked process, nearest neighbors results
# and plot. Also plots pair-correlation function and G_est
hum_spat_fun <- function(data){
  site <- unique(paste0(data$type, data$num))
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
  nn_sd <- sd(nn)
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
  # Plot Gest
  jpeg(paste0("point_process/", site, "_Gest.jpg"))
  # plot(pp, markscale = 0.5,
  #      main = site)
  plot(Gest(pp, correction = "best"),
       main = site,
       xlab = "Distance (m)",
       ylab = "Nearest-neighbour distance CDF")
  dev.off()
  # With envelope
  jpeg(paste0("point_process/", site, "_Gest_envelope.jpg"))
  plot(envelope(pp, Gest, correction = "best", nsim = 19, 
                rank = 1, global = TRUE),
       main = site,
       xlab = "Distance (m)",
       ylab = "Nearest-neighbour distance CDF")
  dev.off()
  
  # Plot average distance between points
  emp <- distmap(pp)
  jpeg(paste0("point_process/", site, "_distmap.jpg"))
  plot(emp, main = site)
  plot(pp, add = TRUE)
  dev.off()
  # Pair correlation = probability of observing a pair of points separated by distance r, divided by corresponding probability for a Poisson process, g(r) = 1 is complete randomness, g(r) > 1 is clustering or attraction, g(r) < 1 is inhibition or regularity
  # Plot pair correlation function
  jpeg(paste0("point_process/", site, "_pcf.jpg"))
  plot(pcf(pp),
       main = site,
       xlab = "Point separation distance (m)",
       ylab = "Nearest neighbor distance probability")
  dev.off()
  
  # # With envelope
  # jpeg(paste0("point_process/", site, "_pcf_envelope.jpg"))
  # plot(envelope(pp, pcf, nsim = 19, rank = 1, global = TRUE),
  #      main = site,
  #      xlab = "Point separation distance (m)")
  # dev.off()
  # Results
  data.frame(nn = nn,
             nn_avg = nn_avg,
             nn_sd = nn_sd,
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

df_hum <- df_n %>%
  mutate(pps = map(data, hum_spat_fun)) %>%
  unnest(pps, .drop = TRUE) %>%
  distinct()
df_hum$pval_text <- ifelse(df_hum$pval_avg < 0.001, 
                                        "<0.001", 
                                        round(df_hum$pval_avg, 3))
# L sites
l <- c("L1", "L2", "L3")

# Predicted normal distributions
normaldens <- df_hum %>%
  group_by(site) %>%
  mutate(
         density = dnorm(nn, nn_avg, nn_sd))
# nn = seq(min(nn), max(nn), length = length(nn)),
# Plot of nn histograms
p_nn <- ggplot(data = dplyr::filter(df_hum,
                                    !(site %in% l)),
               aes(x = nn)) +
  
  geom_histogram(aes(y = ..density..)) +
  geom_line(data = dplyr::filter(normaldens,
                                 !(site %in% l)),
            aes(y = density), colour = "red") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()) + 
  facet_wrap(~site, ncol = 4) +
  geom_text(aes(x = 3.55,
                y = 1.1,
                label = paste("list(mu[NN]==",
                              round(nn_avg, digits = 1),
                              "%+-%",
                              round(nn_sd, digits = 1),
                              "*m",
                              ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(aes(x = 3.55,
                y = 0.95,
                label = paste("list(mu[NN]:mu[exp]==",
                              round(nn_ratio_avg, digits= 2),
                              ")")),
            show.legend = FALSE,
            size = 2,
            parse = TRUE) +
  geom_text(aes(x = 3.7,
                y = 0.8,
                label = paste("p",
                              pval_text)),
            show.legend = FALSE,
            size = 2) +
  ylab("Density") + 
  xlab("Nearest neighbor distance (m)")

ggsave(plot = p_nn,
       filename = "nearest_neighbors.tiff",
       device = "tiff",
       height = 3,
       width = 4,
       units = "in")


# # 
# Fc <- Fest(d2_ppp)
# Fc
# 
# 
# # Nearest neighbor comparison
# Gc <- Gest(d2_ppp)
# Gc
# par(pty = "s")
# plot(Gc)
# 
# # Pairwise distances and K function
# Gc.k <- Kest(d2_ppp)
# Gc.k
# par(pty = "s")
# plot(Gc.k)
# 
# # 
# L <- Lest(d2_ppp)
# plot(L, main = "L function")
# 
# # Pair correlation = probability of observing a pair of points separated by distance r, divided by corresponding probability for a Poisson process, g(r) = 1 is complete randomness, g(r) > 1 is clustering or attraction, g(r) < 1 is inhibition or regularity
# plot(pcf(d2_ppp))
# 
# # J function, Values J(r) > 1 suggest regularity, and J(r) < 1 suggest clustering
# J <- Jest(d2_ppp)
# plot(J)
# 
# # Plot all functions
# plot(allstats(unmark(d2_ppp)))
# 
# # Envelopes
# E <- envelope(d2_ppp, Gest, nsim = 19, rank = 1, global = TRUE)
# E
# plot(E)

# Hummock distribution analysis -------------------------------------------
# Get data in descending rank order
df_h <- df %>%
  dplyr::select(site, area_poly,
                vol, perim_poly) %>%
  dplyr::filter(area_poly > 0.05) %>%
  gather(key = "measure.type", value = "measure",
         -site) %>%
  group_by(site, measure.type) %>%
  mutate(rank = cume_dist(desc(measure)),
         type = str_sub(site, 1, 1))

hist(exp(data))

# exponential model fit with site as dummy variable
e_mod_site <- lm(rank ~ log(measure) * site, 
                 data = filter(df_h,
                               measure.type == "area_poly"))
summary(e_mod_site)
e_mod_all <- lm(log(rank) ~ measure, 
                data = filter(df_h,
                              measure.type == "area_poly"))
summary(e_mod_all)

# Likelihood ratio test for including site in the model
anova(e_mod_site, e_mod_all, test="LRT")


# Model for exponential estimate model at each site
e_mod_site<- function(data) {
  lm(log(rank) ~ measure * site, data = data)
}

# Just group by measure type
mod_site <- df_h %>%
  group_by(measure.type) %>%
  nest() %>%
  mutate(model = purrr::map(data, e_mod_site),
         glance_e = purrr::map(model, glance),
         pval = purrr::map_dbl(glance_e, "p.value"),
         r2 = purrr::map_dbl(glance_e, "adj.r.squared"),
         tidied = purrr::map(model, tidy)) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)


# Definitely should include site in the model (p<<<0.001, adj.r2 = 0.44)
summary_e_mod_site <- tidy(e_mod_site)
write.csv(summary_e_mod_site, "overall_exponential_model_summary.csv")

# Model for exponential estimate model at each site
e_mod <- function(data) {
  lm(log(rank) ~ log(measure), data = data)
}

# Fit exponential line to area, perim, vol, need to filter below 1 a bit
mod_power <- df_h %>%
  group_by(site, measure.type, type) %>%
  filter(measure > 0) %>%
  # dplyr::filter(rank < 0.95) %>%
  nest() %>%
  mutate(model = purrr::map(data, e_mod),
         glance_e = purrr::map(model, glance),
         pval = purrr::map_dbl(glance_e, "p.value"),
         r2 = purrr::map_dbl(glance_e, "adj.r.squared"),
         tidied = purrr::map(model, tidy)) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

write.csv(mod, "exponential_fits_rank.csv")

# Get average fits across sites
mod_sum <- df_h %>%
  dplyr::filter(rank < 0.90,
                type != "B") %>%
  group_by(measure.type) %>%
  nest() %>%
  mutate(model = map(data, e_mod),
         glance_e = map(model, glance),
         pval = map_dbl(glance_e, "p.value"),
         tidied = map(model, tidy)) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate) %>%
  rename(int = `(Intercept)`,
         coef = measure) %>%
  mutate(x = c(10^-1, 10^-3, 10^-4.7),
         y = 0.4)

# refactor things for plotting
df_h$measure.type <- factor(df_h$measure.type,
                            levels = c("perim_poly",
                                       # "area",
                                       "area_poly",
                                       "vol"))
levels(df_h$measure.type) <- c("Perimeter~(m)",
                               # "Area~(m^{2})",
                               "Area~(m^{2})",
                               "Volume~(m^{3})")
# # don't plot these sites
bad <- c("B1", "B3", "B6")

# Plot
p_rank <- ggplot(data = filter(df_h,
                               !(site %in% bad)),
                 aes(x = measure,
                     y = rank)
                ) +
  geom_point(aes(
                 color = site,
                 shape = site),
             alpha = 0.7) +
  scale_shape_manual(name = "Site",
                     labels = c("D1", "D2", "D3", "D4",
                                "T1", "T2", "T3"),
                     values = c(1, 1, 1, 1,
                                16, 16, 16)) +
  scale_color_viridis(name = "Site",
                      labels = c("D1", "D2", "D3", "D4",
                                 "T1", "T2", "T3"),
                      discrete = TRUE) +
  stat_smooth(method = "lm",
              formula = y ~ exp(x),
              se = FALSE,
              show.legend = FALSE) +
  geom_text(data = mod_sum,
            aes(x = x,
                y = y,
                label = paste("list(",
                              round(exp(int), digits = 2), 
                              "*e^{",
                              round(coef, digits = 2),
                              "*x",
                              "}",
                              ")")),
            show.legend = FALSE,
            size = 2.5,
            parse = TRUE) +
  theme_bw() +
  xscale("log10", .format = TRUE) +
  scale_y_log10(
         labels = trans_format('log10', math_format(10 ^ .x)),
         limits = c(10 ^ -3.2, 10 ^ 0.2),
         breaks = c(10^-3, 10 ^ -2, 10 ^ -1, 10 ^ 0)
         ) +
  annotation_logticks(base = 10, sides = "tlbr") + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1,1,1,1), "cm"),
    aspect.ratio = 1,
    legend.position = c(0.08, 0.39),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.background = element_rect(
      fill = "gray90",
      linetype = "solid",
      colour = "black")) + 
  xlab("") +
  ylab(expression("P (" * X>=x * ")")) + 
  facet_wrap(~measure.type, ncol = 3, scales = "free_x", 
             labeller = label_parsed)
ggsave(p_rank, filename = "Hummock_properties_rank.tiff",
       device = "tiff",
       dpi = 300,
       width = 6,
       height = 4)



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



# Power law analysis ------------------------------------------------------
xmins <- df_h %>%
  group_by(site, measure.type) %>%
  distinct(measure) %>%
  arrange(desc(measure))
data <- df_h %>%
  filter(site == "D1",
         measure.type == "area_poly") %>%
  pull(measure)
xmins <- df_h %>%
  filter(site == "D1",
         measure.type == "area_poly") %>%
  distinct(measure) %>%
  pull(measure)

m_bl_ln = conlnorm$new(data)
m_bl_ln$setPars(estimate_pars(m_bl_ln))
est = estimate_xmin(m_bl_ln)
m_bl_ln$setXmin(est)
lines(m_bl_ln, col=3, lwd=2)

m_bl_exp = conexp$new(data)
m_bl_exp$setPars(estimate_pars(m_bl_exp))
est2 = estimate_xmin(m_bl_exp)
m_bl_exp$setXmin(est)
lines(m_bl_exp, col=4, lwd=2)

com <- compare_distributions(m_bl_exp, m_bl_ln)
com$p_two_sided



m1 = conpl$new(data)
m1$setXmin(estimate_xmin(m1))
m2 = conlnorm$new(data)
m2$setXmin(m1$getXmin())
m2$setPars(estimate_pars(m2))
plot(m2, ylab="CDF")
lines(m1)
lines(m2, col=2, lty=2)
comp = compare_distributions(m1, m2)
comp$p_two_sided






dat = numeric(length(xmins))
z = sort(data)

for (i in 1:length(xmins)){
  xmin = xmins[i] # choose next xmin candidate
  z1 = z[z>=xmin] # truncate data below this xmin value
  n = length(z1)
  a = 1+ n*(sum(log(z1/xmin)))^-1 # estimate alpha using direct MLE
  cx = (n:1)/n # construct the empirical CDF
  cf = (z1/xmin)^(-a+1) # construct the fitted theoretical CDF
  dat[i] = max(abs(cf-cx)) # compute the KS statistic
}
D = min(dat[dat>0],na.rm=TRUE) # find smallest D value
xmin = xmins[which(dat==D)] # find corresponding xmin value
z = data[data>=xmin]
z = sort(z)
n = length(z)
alpha = 1 + n*(sum(log(z/xmin)))^-1 # get corresponding alpha estimate
library(gsl)
library(numDeriv)

dpowerlaw <- function(x, alpha=2, xmin=1, log=F) {
  if (log)
    log(alpha-1) - log(xmin) - alpha * log(x / xmin)
  else
    ((alpha - 1) / xmin) * ((x / xmin) ^ (-alpha))
}
ppowerlaw <- function(q, alpha=2, xmin=1, lower.tail=T, log.p = F) {
  p <- (q / xmin) ^ (- alpha + 1)
  if (lower.tail)
    p <- 1-p
  if (log.p)
    p <- log(p)
  p
}
qpowerlaw <- function(p, alpha=2, xmin=1, lower.tail=T, log.p = F) {
  if (!lower.tail)
    p <- 1-p
  if (log.p)