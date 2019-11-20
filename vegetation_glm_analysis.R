# 
# Author: Jake Diamond
# Purpose: To analyze vegetation field data categorically
# Date: November 16, 2016
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(broom)
library(Hmisc)
library(gridExtra)
library(lme4)
library(MASS)
library(tidyverse)

# ggplot theme
theme_set(theme_bw() + 
            theme(panel.grid = element_blank()))

# Load veg data
df <- read_rds("veg_data.rds") %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              ifelse(site == "B6",
                                     "L3",
                                     site))))

# Get data into format for vegan diversity function, and calc richness
df_r <- df %>%
  mutate(tot = ifelse(is.na(number),
                            percent,
                            number)) %>%
  group_by(site, point, plot, depth, hu.ho, moss, notree) %>%
  dplyr::summarize(rich = sum(tot > 0, na.rm = T))

# Get elevation data
elev <- read_csv("relative_elevations_all_v6.csv") %>%
  dplyr::filter(!is.na(point)) %>%
  dplyr::select(-X1, -(id:PA_poly))
# Use the z coord at each site's well as our datum
elev <- elev %>%
  dplyr::filter(point == "well") %>%
  transmute(site = site,
            z_well = z) %>%
  right_join(elev) %>%
  mutate(z = z - z_well,
         zrel = z - mean)

# Join data
dfr <- elev %>%
  dplyr::select(site, point, z, zrel) %>%
  right_join(df_r, by = c("site", "point")) %>%
  mutate(tree = ifelse(hu.ho == "hummock" &
                         is.na(notree),
                       1,
                       0))
dfr$depth <- as.numeric(dfr$depth)

# summary plots
rich <- ggplot(data = dfr, 
               aes(x = site, y = rich))+
  stat_summary(fun.data = mean_cl_boot) + 
  xlab("") + ylab("Sample richness")
rich
pts <- ggplot(data = dfr, 
              aes(x = zrel,y = rich,
                  colour = site, group = site))+
  geom_point(alpha = 0.2) + 
  scale_color_viridis_d(name = "Site") + 
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  xlab("Relative elevation (m)") +
  ylab("Sample richness")
pts
ggsave(rich,
       filename = "veg_models/richness_summary.tiff",
       device = "tiff",
       dpi = 300)
ggsave(pts,
       filename ="veg_models/richness_by_elev.tiff",
       device = "tiff",
       dpi = 300)

# Basic model
mod0 <- glmer(rich ~ z + (1|site),
              family = "poisson",
              data = dfr,
              control = glmerControl(optimizer="bobyqa",
                                     check.conv.grad = .makeCC("warning",
                                                               2e-3)))
mod0
# Random slope effect
modrs <- glmer(rich ~ z + (1+z|site),
              family = "poisson",
              data = dfr,
              control = glmerControl(optimizer="bobyqa",
                                     check.conv.grad = .makeCC("warning",
                                                               2e-3)))
# Compare those two
anova(mod0, modrs)

# Random slope is better(p = 0.0078)
# Include plot nesting
modrsp <- glmer(rich ~ z + (1 + z|site/plot),
                     family = "poisson",
                     data = dfr,
                     control = glmerControl(optimizer="bobyqa",
                                          check.conv.grad = .makeCC("warning",
                                                                  2e-3)))
anova(modrs, modrsp)
# Not going to include plot (p = 0.04)
# Include moss difference
modrsm <- glmer(rich ~ z + moss + (1 + z|site),
                     family = "poisson",
                     data = dfr,
                     control = glmerControl(optimizer="bobyqa",
                                            check.conv.grad = .makeCC("warning",
                                                                      2e-3)))
anova(modrsm, modrs)
# Definitely include moss (p<<<0.0001)
summary(modrsm)

# Diagnostics
png(paste0("veg_models/diagnostics.png"),
    width = 800, height = 600)
grid.arrange(plot(mod0,type=c("p","smooth")),
                  plot(mod0,sqrt(abs(resid(.)))~fitted(.),
                       type=c("p","smooth"),
                       ylab=expression(sqrt(abs(resid)))))
dev.off()

# Random effects
png(paste0("veg_models/ranef.png"),
    width = 800, height = 600)
dotplot(ranef(modrsm,condVar=TRUE),
              lattice.options=list(layout=c(1,2)))
dev.off()

# Prediction without random effects
pframe <- expand.grid(
                      moss = c(0, 1),
                      z = seq(min(dfr$z, na.rm = T), 
                              max(dfr$z, na.rm = T), 0.05))
pred1 <- predict(modrsm, re.form = NA,
                 newdata = pframe,
                 type = "response")
# Bootstrap for confidence intervals
set.seed(101)
bb <- bootMer(modrsm,
              FUN=function(x)
                predict(x, re.form = NA,
                        newdata = pframe,
                        type = "response"),
              nsim = 100)
predboot1.CI <- t(apply(bb$t, 2,
                        quantile,
                        c(0.025,0.975),
                        na.rm = TRUE))

pframe2 <- cbind(pframe,
                 rich = pred1,
                 setNames(as.data.frame(predboot1.CI),
                          c("lwr_boot", "upr_boot")))

# Plot of predictions
predp <- ggplot(data = pframe2,
                aes(x = z,
                    y = rich,
                    linetype = as.factor(moss))) +
  geom_line() +
  geom_ribbon(aes(ymin=lwr_boot,ymax=upr_boot),
              alpha=0.05) +
  # geom_point(data = dfr,
  #            aes(color = site)) + 
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  scale_linetype_manual(name = "",
                        breaks = c(0, 1),
                        values = c("solid", "dashed"),
                        labels = c("Vascular", "Moss")) + 
  ylab(expression("Predicted richness (species 0.25"*m^{-2}*")")) + 
  xlab("Relative elevation above water table (m)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.8))
predp

ggsave(predp, 
       filename = "veg_models/veg_predictions_ribbon.tiff",
       device = "tiff",
       dpi = 300,
       width = 4,
       height = 3,
       units = "in")

# Tidied model data
tidied <- tidy(modrsm)
ranef(modrsm)
write.csv(tidied, "veg_models/veg_glm_summary.csv")
