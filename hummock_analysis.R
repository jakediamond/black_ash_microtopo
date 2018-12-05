# 
# Author: Jake Diamond
# Purpose: Analysis of hummocks, point processes
# Date: December 5, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(spatstat)
library(tidyverse)
library(broom)
library(MASS)
library(scales)

# Get hummock data (x, y, area, vol, perim)
df <- read.csv("Lidar/site_hum_stats.csv")
df_d2 <- df[df$site == "T1", ]
df_d2.big <- dplyr::filter(df_d2, area > 0.3)
d2_ppp <- ppp(x = df_d2$x,
              y = df_d2$y,
              xrange = c(-15, 15),
              yrange = c(-20, 10),
              marks = df_d2$area)
unitname(d2_ppp) <- c("meter", "meter")

# Plot marked spatial point process
plot(d2_ppp)

# Summarize data and get the estimated intensity
summary(d2_ppp)
lamb <- summary(d2_ppp)$intensity
hist(marks(d2_ppp))

# Quadrat counting to test for inhomogeneity of intensity
Q <- quadratcount(d2_ppp,
                  nx = 4,
                  ny = 6)

# Plot quadrats
plot(d2_ppp)
plot(Q, add = TRUE, cex = 2)

# Kernel density (for intensity)
den <- density(d2_ppp, sigma = 10)
plot(den)
persp(den)
contour(den)

# Create a window
W <- owin(xrange = c(0, 2),
          yrange = c(0, 3))
d2_ppp[W]

# Quadrat test chi-square for complete spatial randomness
M <- quadrat.test(d2_ppp, 
                  nx = 5,
                  ny = 5)
plot(d2_ppp)
plot(M, add = TRUE)
M
# KS test for CSR
KS <- cdf.test(d2_ppp, "x")
plot(KS)
KS
# Test against Poisson homogeneity
lambda <- function(x, y){
  100 * (x + y)
}
plot(rpoispp(lambda))

# Homogenous poisson model
ppm(unmark(d2_ppp), ~1)

# Inhomogenous poisson model with intensity that is log-linear in the cartestian coordinates
fit <- ppm(unmark(d2_ppp), ~x + y)
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








df_d2.h <- df_d2 %>%
  gather(key = "measure.type", value = "measure",
         -x,
         -y,
         -site) %>%
  group_by(measure.type) %>%
  mutate(rank = cume_dist(desc(measure)))
 

# Fit exponential line
mod <- tidy(lm(rank ~ exp(area), data = df_d2.h))

# Plot
p_exp2 <- ggplot(data = df_d2.h,
                aes(x = measure,
                    y = rank)) +
  geom_point(size = 3,
             shape = 1) +
  stat_smooth(formula = y ~ exp(x)) +
  facet_wrap(~measure.type) +
  # annotate("text",
  #          x = 0.001,
  #          y = 0.1,
  #          label = "p = 1.8E-18",
  #          size = 6) + 
  theme_bw() +
  scale_x_log10(
    labels = trans_format('log10', math_format(10 ^ .x))
    # ,
    # limits = c(10 ^ -3.3, 10 ^ -0.3),
    # breaks = c(10 ^ -3, 10 ^ -2, 10 ^ -1)
    ) +
  scale_y_log10(
    labels = trans_format('log10', math_format(10 ^ .x))
    # ,
    # limits = c(10 ^ -2, 10 ^ 0),
    # breaks = c(10 ^ -2, 10 ^ -1, 10 ^ 0)
    ) +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 18, colour = "black"),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.12, 0.17),
    plot.title = element_text(lineheight = .8, face = "bold"),
    legend.text = element_text(size = 16),
    legend.background = element_rect(
      fill = "gray90",
      size = 1,
      linetype = "solid",
      colour = "black"))
  # xlab(expression("Hummock volume, v (" * m ^ 3 * ")")) +
  # ylab(expression("P (" * V>=v* ")"))
p_exp2
p_exp

ggsave(p_exp, filename = "D2_hummock_exponential.tiff",
       device = "tiff",
       dpi = 600)
