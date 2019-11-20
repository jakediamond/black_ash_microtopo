# 
# Author: Jake Diamond
# Purpose: Analysis of hummocks, point processes
# Date: December 5, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")
setwd("C:/Users/jake.diamond/Dropbox/Projects/EAB/Data")

# Load Libraries
library(raster)
library(MASS)
library(spatstat)
library(gridExtra)
library(ncf)
library(nlme) 
library(Hmisc)
library(tidyverse)
library(broom)
library(scales)
library(sf)
library(lme4)
library(poweRlaw)
library(fitdistrplus)
library(ggpubr)
library(furrr)

# Get hummock data (x, y, area, vol, perim), rename sites
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

# Get summary of all hummock data
df_h <- df %>%
  dplyr::select(site, area_poly,
                vol, perim_poly, PA_poly) %>%
  dplyr::filter(between(area_poly, 0.05, 100)) %>%
  dplyr::gather(key = "measure.type", value = "measure",
         -site) %>%
  filter(measure > 0)

df_h_sum <- df_h %>%
  group_by(site, measure.type) %>%
  mutate(ln = log(measure)) %>%
  dplyr::summarize(mean = exp(mean(ln) + sd(ln)^2 / 2),
            sd = sqrt((exp(sd(ln)^2) -1) * exp(2 * mean + sd(ln)^2)))


# Quickly plot hummock perimeter vs area
pa_p <- ggplot(data = filter(df,
                             between(area_poly, 0.05, 30)),
               aes(x = perim_poly,
                   y = area_poly)) +
  geom_point(aes(color = site)) +
  stat_smooth(method = "glm",
              formula = y ~ x + I(x^2)) +
  scale_color_viridis_d() +
  theme_bw()
pa_p

# Estimate how much additional surface area is provided by 
# hummocks. Using cone geometry estimate
df_sa <- df %>%
  filter(area > min_hum) %>%
  transmute(site = site,
            site_area = site_area,
            sa = zmeann * perim_poly,
            sa_cone = 3.14 * sqrt(3*vol*zmaxn/3.14+(3*vol/(3.14*zmaxn))^2)) %>%
  group_by(site) %>%
  dplyr::summarize(site_area = mean(site_area),
            sa = sum(sa),
            sa_cone = sum(sa_cone),
            sa_add = sa / site_area,
            sa_cone_add = sa_cone / site_area)

# Estimate how much of estimated wetland volume is taken 
# up by hummocks. Wetland volume is estimated as the mean 80th 
# quantile of hummock height in a site * wetland area
df_vol <- df %>%
  group_by(site) %>%
  dplyr::summarize(site_ht = mean(z80n),
            site_vol = site_ht * mean(site_area),
            hum_vol = sum(vol),
            hum_ratio = hum_vol / site_vol)

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
  # Plot Gest, not doing this anymore because already done, not in paper
  # jpeg(paste0("point_process/", site, "_Gest.jpg"))
  # # plot(pp, markscale = 0.5,
  # #      main = site)
  # plot(Gest(pp, correction = "best"),
  #      main = site,
  #      xlab = "Distance (m)",
  #      ylab = "Nearest-neighbour distance CDF")
  # dev.off()
  # # With envelope
  # jpeg(paste0("point_process/", site, "_Gest_envelope.jpg"))
  # plot(envelope(pp, Gest, correction = "best", nsim = 19, 
  #               rank = 1, global = TRUE),
  #      main = site,
  #      xlab = "Distance (m)",
  #      ylab = "Nearest-neighbour distance CDF")
  # dev.off()
  # 
  # # Plot average distance between points
  # emp <- distmap(pp)
  # jpeg(paste0("point_process/", site, "_distmap.jpg"))
  # plot(emp, main = site)
  # plot(pp, add = TRUE)
  # dev.off()
  # # Pair correlation = probability of observing a pair of points separated by distance r, divided by corresponding probability for a Poisson process, g(r) = 1 is complete randomness, g(r) > 1 is clustering or attraction, g(r) < 1 is inhibition or regularity
  # # Plot pair correlation function
  # jpeg(paste0("point_process/", site, "_pcf.jpg"))
  # plot(pcf(pp),
  #      main = site,
  #      xlab = "Point separation distance (m)",
  #      ylab = "Nearest neighbor distance probability")
  # dev.off()
  
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
  mutate(pps = future_map(data, hum_spat_fun)) %>%
  unnest(pps, .drop = TRUE) %>%
  distinct()
df_hum$pval_text <- ifelse(df_hum$pval_avg < 0.001, 
                                        "<0.001", 
                                        paste0("=", round(df_hum$pval_avg, 3)))
# L sites
l <- c("L1", "L2", "L3")

# Predicted normal distributions
normaldens <- df_hum %>%
  group_by(site) %>%
  mutate(
         density = dnorm(nn, nn_avg, nn_sd))
# nn = seq(min(nn), max(nn), length = length(nn)),
# Plot of nn histograms
p_nn <- ggplot(data = df_hum,
               aes(x = nn)) +
  geom_histogram(aes(y = ..density..)) +
  geom_line(data = normaldens,
            aes(y = density), colour = "red") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        strip.text.x = element_text(margin = margin(0,0,0,0, "cm"),
                                    size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)) + 
  # scale_x_continuous(limits = c(0, 10),
  #                    breaks = seq(0, 10, 2.5)) +
  facet_wrap(~site, ncol = 4) +
  geom_text(aes(x = 6.7,
                y = 1.1,
                label = paste("list(mu[NN]==",
                              round(nn_avg, digits = 1),
                              "%+-%",
                              round(nn_sd, digits = 1),
                              "*m",
                              ")")),
            show.legend = FALSE,
            size = 1.8,
            parse = TRUE) +
  geom_text(aes(x = 7,
                y = 0.95,
                label = paste("list(mu[NN]:mu[exp]==",
                              round(nn_ratio_avg, digits= 2),
                              ")")),
            show.legend = FALSE,
            size = 1.8,
            parse = TRUE) +
  geom_text(aes(x = 7.15,
                y = 0.8,
                label = paste("p",
                              pval_text)),
            show.legend = FALSE,
            size = 1.8) +
  ylab("Density") + 
  xlab("Nearest neighbor distance (m)")
p_nn
ggsave(plot = p_nn,
       filename = "Figures/nearest_neighbors_l_publication2.tiff",
       device = "tiff",
       dpi = 300,
       height = 3.5,
       width = 5,
       units = "in")

# Hummock distribution analysis -------------------------------------------
# Get data in descending rank order
df_h_desc <- df_h %>%
  group_by(site, measure.type) %>%
  mutate(rank = cume_dist(desc(measure)),
         type = str_sub(site, 1, 1))
psd <- filter(df_h_desc,
              measure.type == "area_poly",
              site == "D1") %>%
  pull(measure)
xmin <- 0.1
library(spatialwarnings)
lib
# Fit a model for everyone. Note that we store the results of the variables 
# here so we reuse the previous expo and rate for tpl fitting. 
plfit  <- pl_fit(psd, xmin)
expfit <- exp_fit(psd, xmin)
tplfit <- tpl_fit(psd, xmin)
lnormfit <- lnorm_fit(psd, xmin)
models <- list(pl = plfit, tpl = tplfit, exp = expfit,
               lnorm = lnormfit)

models <- lapply(models, as.data.frame)
models <- do.call(bind_rows, models)
row.names(models) <- models[ ,'type']

# Compute AICs
models[ ,'AIC']  <- get_AIC(models[ ,'ll'],  models[ ,'npars'])
models[ ,'AICc'] <- get_AICc(models[ ,'ll'], models[ ,'npars'], length(psd))
models[ ,'BIC']  <- get_BIC(models[ ,'ll'],  models[ ,'npars'], length(psd))

best_by <- "AIC"
# We need to remove NA's here as sometimes one of the fits fails and its "best"
# column is NA. We do not consider failed fit as good solutions. 
models[ ,'best'] <- models[ ,best_by] == min(models[ ,best_by], na.rm = TRUE)
models[ ,'best'] <- models[ ,'best'] & ! is.na(models[ ,'best']) 

# Add an xmin column 
models[ ,"xmin_fit"] <- xmin

# Reorganize columns 
models <- models[ ,table_names]

get_AIC <- function(ll, k) { 
  2*k - 2*ll 
}

get_AICc <- function(ll, k, n) { 
  2*k - 2*ll + (2*k*(k+1))/(n-k-1)
}

get_BIC <- function(ll, k, n) { 
  2*k*log(n) - 2*ll 
}
newdata <- unique( round(10^(seq(0, log10(max(psd)), 
                                 length.out = 200))) )
vals_pred <- data.frame()
lnorm = pdislnorm(newdata, 
                  models["lnorm", "meanlog"], 
                  models["lnorm", "sdlog"],  
                  models["lnorm", "xmin_fit"])
vals_pred <- rbind(vals_pred, 
                   data.frame(type = "lnorm", patchsize = newdata, 
                              y = lnorm))
vals_pred <- vals_pred[ vals_pred[ ,'y'] >= min(psd), ] 
ggplot() + 
  # scale_y_log10() +
  # scale_x_log10() + 
  xlab('Patch size') + 
  ylab('Frequency (P>=x)') + 
  geom_point(aes(x = measure, y = rank), data = filter(df_h_desc,
                                                          measure.type == "area_poly",
                                                          site == "D1"))  +
  stat_ecdf(color="red")
glance(fitdistr(data, "lognormal"))
glance(fitdistr(data, "exponential"))
glance(fitdistr(data, "gamma"))
# Nest data
dfn <- df_h_desc %>%
  mutate(mt = measure.type) %>%
  group_by(site, mt) %>%
  nest()

test <- dfn %>%
  mutate(ln = map(data, ~fitdistr(.$measure, "lognormal")),
         e = map(data, ~fitdistr(.$measure, "exponential")),
         g = map(data, ~fitdistr(.$measure, "gamma")),
         lna = map(ln, glance),
         ea = map(e, glance),
         ga = map(g, glance)) %>%
  unnest(lna, ea, ga) %>%
  gather(AIC, value, starts_with("AIC")) %>%
  group_by(site, mt) %>%
  filter(value == min(value)) %>%
  mutate(ml = map(ln, pluck, 1, 1),
         sdl = map(ln, pluck, 1, 2)) %>%
  dplyr::select(site, mt, ml, sdl) %>%
  unnest()

ggplot() + 
  # scale_y_log10() +
  # scale_x_log10() +
  xlab('Patch size') + 
  ylab('Frequency (P>=x)') + 
  # geom_point(aes(x = measure, y = rank), data = filter(df_h_desc,
  #                                                      measure.type == "area_poly",
  #                                                      site == "D1"))  +
  stat_function(fun = plnorm, 
                args = list(q = 1:20, meanlog = -0.5894042
, sdlog = 1.3636570,
                            lower.tail  =FALSE),
                  color="red")

# Function that does lmer analysis
lmfun <- function(data){
  theme_set(theme_bw())
  mt <- unique(data$measure.type)
  # summary plots
  sitecon <- ggplot(data,aes(x=site,y=measure))+
    stat_summary(fun.data=mean_cl_boot) + 
    xlab("") + ylab("Value") + 
    ggtitle(mt)
  ggsave(sitecon,
         filename = paste0("hummock_ranks/",
                           mt, 
                           "_summary.tiff"),
         device = "tiff",
         dpi = 300)

  # initial model
  mod0 <- lme(rank ~ log(measure),
              data = data,
              method = "REML",
              random = ~ 1 + log(measure)|site,
              control=list(maxIter=10000, niterEM=10000))
  
  # Diagnostics
  png(paste0("hummock_ranks/", mt, "_diagnostics.png"),
      width = 800, height = 600)
  d <- grid.arrange(plot(mod0,type=c("p","smooth")),
                    plot(mod0,sqrt(abs(resid(.)))~fitted(.),
                         type=c("p","smooth"),
                         ylab=expression(sqrt(abs(resid)))),
                    plot(mod0,resid(.,type="pearson")~log(measure),
                         type=c("p","smooth")),
                    qqnorm(mod0,abline=c(0,1)))
  print(d)
  dev.off()
  mod0
}

lm4fun <- function(data){
  mt <- unique(data$measure.type)
  data$measure
  mod0r <- glmer(rank ~ log(measure) + (log(measure)|site),
                REML = TRUE,
                data = data,
                family = gaussian(link = "log"))
  # Random effects
  png(paste0("hummock_ranks/", mt, "_ranef.png"),
      width = 800, height = 600)
  re <- dotplot(ranef(mod0r,condVar=TRUE),
                lattice.options=list(layout=c(1,2)))
  print(re)
  dev.off()
  mod0r
}

predfun <- function(data, mod){
  mt <- unique(data$measure.type)
  # Prediction
  pframe <- expand.grid(site = as.factor(data$site),
                        measure = seq(from = min(data$measure), 
                                      to = max(data$measure),
                                      length.out = 25))
  pred1 <- predict(mod,re.form=NA,
                   newdata=pframe,
                   type="response")
  set.seed(101)
  bb <- bootMer(mod,
                FUN=function(x)
                  predict(x,re.form=NA,
                          newdata=pframe,
                          type="response"),
                nsim=500)
  predboot1.CI <- t(apply(bb$t, 2,
                          quantile,
                          c(0.025,0.975),
                          na.rm=TRUE))
  
  pframe2 <- cbind(pframe,rank=pred1,
                   setNames(as.data.frame(predboot1.CI),
                            c("lwr_boot","upr_boot")))
  
}

# Apply linear mixed effect model functions to data
lmout <- dfn %>%
  group_by(mt) %>%
  mutate(
    lmres = map(data, lmfun),
    lm4res = map(data, lm4fun)
    ,preds = map2(data, lm4res, predfun)
  )
# Tidied model data
tidied <- lmout %>%
  dplyr::select(-data) %>%
  mutate(tid = map(lm4res, tidy)
         ,ano = map(lmres, anova),
         pval = unlist(map(ano, pluck, 4, 2))
  ) %>%
  dplyr::select(-lm4res, -preds, -lmres, -ano) %>%
  unnest()

write.csv(tidied, "hummock_ranks/cdf_ranef_models_lognormal.csv")


# predictions with confidence intervals (no random effects)
preds <- lmout %>%
  dplyr::select(-data, -lm4res, -lmres) %>%
  unnest() %>%
  dplyr::select(-site) %>%
  distinct()


write.csv(preds, "predictions_ranef_rank.csv")
# Plot of predictions
ggplot(data = preds,
                aes(x = measure,
                    y = rank)) +
  facet_wrap(~mt, scale = "free_y")+
  geom_line() +
  geom_ribbon(aes(ymin=lwr_boot,ymax=upr_boot),
              alpha=0.05) +
  # xscale("log10", .format = TRUE) +
  yscale("log10", .format = TRUE) +
  scale_color_viridis_d() + 
  ylab(expression("Predicted concentration (mg"~L^{-1}*"or -)")) +
  xlab("Relative elevation (m)")
  



# Get average fits across sites
mod_sum <- df_h_desc %>%
  dplyr::filter(type != "L") %>%
  group_by(measure.type) %>%
  nest() %>%
  mutate(model = map(data, ln_mod),
         glance_e = map(model, glance),
         r2 = purrr::map_dbl(glance_e, "adj.r.squared"),
         pval = map_dbl(glance_e, "p.value"),
         tidied = map(model, tidy)) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  filter(measure.type != "PA_poly") %>%
  spread(term, estimate) %>%
  rename(int = `(Intercept)`,
         coef = `log(measure)`) %>%
  mutate(x = c(10^-1, 10^0, 10^-2.5),
         y = 0.3)
# Only plot area_poly
df_h2 <- filter(df_h_desc, 
                !(measure.type %in% c("PA_poly", "area")))
# refactor things for plotting
df_h2$measure.type <- factor(df_h2$measure.type,
                            levels = c("perim_poly",
                                       # "area",
                                       "area_poly",
                                       "vol"))
levels(df_h2$measure.type) <- c("Perimeter~(m)",
                               # "Area~(m^{2})",
                               "Area~(m^{2})",
                               "Volume~(m^{3})")

mod_sum$measure.type <- factor(mod_sum$measure.type,
                             levels = c("perim_poly",
                                        # "area",
                                        "area_poly",
                                        "vol"))
levels(mod_sum$measure.type) <- c("Perimeter~(m)",
                                # "Area~(m^{2})",
                                "Area~(m^{2})",
                                "Volume~(m^{3})")
# Plot
p_rank <- ggplot(data = df_h_desc,
                 aes(x = measure,
                     y = rank)
                ) +
  geom_point(aes(
                 color = site),
             alpha = 0.7) +
  # scale_shape_manual(name = "Site",
  #                    labels = c("D1", "D2", "D3", "D4",
  #                               "T1", "T2", "T3"),
  #                    values = c(1, 1, 1, 1,
  #                               16, 16, 16)) +
  scale_color_viridis(name = "Site",
                      discrete = TRUE) +
  stat_smooth(method = "glm",
              formula = y ~ exp(x),
              se = FALSE,
              show.legend = FALSE) +
  # geom_text(data = mod_sum,
  #           aes(x = x,
  #               y = y,
  #               label = paste("list(y==",
  #                             round(coef, digits = 2),
  #                             "*ln(x)+",
  #                             round(int, digits = 2),
  #                             ")")),
  #           show.legend = FALSE,
  #           size = 2,
  #           parse = TRUE) +
  # geom_text(data = mod_sum,
  #           aes(x = x,
  #               y = y - 0.1,
  #               label = paste("list(R^2==",
  #                             round(r2, digits = 2),
  #                             ")")),
  #           show.legend = FALSE,
  #           size = 2,
  #           parse = TRUE) +
  theme_bw() +
  xscale("log10", .format = TRUE) +
  yscale("log10", .format = TRUE) +
  # scale_y_log10(
  #        labels = trans_format('log10', math_format(10 ^ .x)),
  #        limits = c(10 ^ -3.2, 10 ^ 0.2),
  #        breaks = c(10^-3, 10 ^ -2, 10 ^ -1, 10 ^ 0)
  #        ) +
  # annotation_logticks(base = 10, sides = "tlbr") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(margin = margin(0,0,0,0, "cm")),
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
  guides(color = guide_legend(ncol = 2)) +
  # facet_wrap(~measure.type, ncol = 3, scales = "free_x", 
  #            labeller = label_parsed)
  facet_grid(type~measure.type, scales = "free_x")
p_rank

ggsave(p_rank, filename = "Hummock_properties_rank_ln.tiff",
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
hist(data)
fit.exp <- fitdist(data, "exp")
plot(fit.exp)
fit.exp <- fitdist(data, "lnorm")
plot(fit.exp)

m_bl_ln = conlnorm$new(data)
m_bl_ln$setPars(estimate_pars(m_bl_ln))
est = estimate_xmin(m_bl_ln)
m_bl_ln$setXmin(est)
plot(m_bl_ln)
lines(m_bl_ln, col=3, lwd=2)

m_bl_exp = conexp$new(data)
m_bl_exp$setPars(estimate_pars(m_bl_exp))
est2 = estimate_xmin(m_bl_exp)
m_bl_exp$setXmin(est2)
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