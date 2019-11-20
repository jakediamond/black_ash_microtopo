# 
# Author: Jake Diamond
# Purpose: To analyze soil chem data with linear mixed effect
# models of elevation
# Date: February 25, 2019
# 

# Set Working Directory
setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(broom)
library(vegan)
library(gridExtra)
library(lme4)
library(MASS)        
library(nlme) 
library(tidyverse)
# ggplot theme
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(),
                  strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))))
# Load chem data
df <- read_rds("soil_chem_long")

# Some outlier functions
# Z score
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}
# Median absolute deviation
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}
# Tukey's fences
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

# Apply outlier functions, get rid of outliers, and clean
df <- df %>%
  group_by(site, solute) %>%
  mutate(insz = isnt_out_z(conc),
         insmad = isnt_out_mad(conc),
         instuk = isnt_out_tukey(conc),
         outlier = ifelse((insz + insmad +instuk) < 2,
                          1,
                          0)) %>%
  dplyr::filter(outlier == 0) %>%
  mutate(tree = ifelse(hu.ho == "hummock" &
                         notree == 0,
                       1,
                       0),
         zrel = z - z_well,
         zw = zrel - mean)

# Nest data
dfn <- df %>%
  group_by(site, solute) %>%
  mutate(conc_s = scale(conc)) %>%
  ungroup() %>%
  mutate(sol = solute) %>%
  group_by(solute) %>%
  nest()

# Function that does lmer analysis
lmfun <- function(data){
  theme_set(theme_bw())
  sol <- unique(data$sol)
  # summary plots
  sitecon <- ggplot(data,aes(x=site,y=conc_s))+
    stat_summary(fun.data=mean_cl_boot) + 
    xlab("") + ylab("Scaled concentration (mg/L)")
  conpts <- ggplot(data, 
                   aes(x=zw,y=conc_s,colour=site,group=site))+
    geom_point()+scale_color_viridis_d() + 
    xlab("Relative elevation (m)") +
    ylab("Scaled concentration (mg/L)")
  ggsave(sitecon,
         filename = paste0("soils_models/",
                           sol, 
                           "_scaled_concentration_summary.tiff"),
         device = "tiff",
         dpi = 300)
  ggsave(conpts,
         filename = paste0("soils_models/",
                           sol, 
                           "_scaled_concentrations_by_elev.tiff"),
         device = "tiff",
         dpi = 300)
  # initial model
  mod0 <- lme(conc_s ~ zw,
              data = data,
              method = "REML",
              random = list(site = pdDiag(~ zw)),
              control=list(maxIter=10000, niterEM=10000))

  # Diagnostics
  png(paste0("soils_models/", sol, "scaled_diagnostics.png"),
      width = 800, height = 600)
  d <- grid.arrange(plot(mod0,type=c("p","smooth")),
               plot(mod0,sqrt(abs(resid(.)))~fitted(.),
                    type=c("p","smooth"),
                    ylab=expression(sqrt(abs(resid)))),
               plot(mod0,resid(.,type="pearson")~z,
                    type=c("p","smooth")),
               qqnorm(mod0,abline=c(0,1)))
  print(d)
  dev.off()
 mod0
}

lm4fun <- function(data){
  sol <- unique(data$sol)
  mod0r <- lmer(conc_s ~ zw + (zw||site),
                REML = TRUE,
                data = data
                )
  # Random effects
  png(paste0("soils_models/", sol, "scaled_ranef.png"),
      width = 800, height = 600)
  re <- dotplot(ranef(mod0r,condVar=TRUE),
          lattice.options=list(layout=c(1,2)))
  print(re)
  dev.off()
  mod0r
}

predfun <- function(data, mod){
  sol <- unique(data$sol)
  # Prediction
  pframe <- expand.grid(site = as.factor(data$site),
                        zw = seq(min(data$zw), max(data$zw), 0.1))
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

  pframe2 <- cbind(pframe,conc_s=pred1,
                   setNames(as.data.frame(predboot1.CI),
                            c("lwr_boot","upr_boot")))

}

# Apply linear mixed effect model functions to data
lmout <- dfn %>%
  group_by(solute) %>%
  mutate(
    lmres = map(data, lmfun)
    ,lm4res = map(data, lm4fun)
    ,preds = map2(data, lm4res, predfun)
    ,re = map(lmres, ranef)
    ,fe = map(lmres, fixef)
  )
map(lmout$re, plot)

# Get ran effects
es <- lmout %>%
  dplyr::select(solute, re) %>%
  mutate(ran = map(re, pluck, 1)) %>%
  unnest(ran) %>%
  mutate(site = unique(df$site)) %>%
  spread(solute, ran)
write.csv(es, "soils_models/random_effects_lm3_int.csv")



# get fixed effects
fs <- lmout %>%
  dplyr::select(solute, fe) %>%
  unnest(fe) %>%
  mutate(type = c("int", "slope")) %>%
  spread(type, fe)
write.csv(fs, "soils_models/fixed_effects.csv")

# predictions with confidence intervals (no random effects)
preds <- lmout %>%
  dplyr::select(-data, -lm4res, -lmres) %>%
  unnest() %>%
  dplyr::select(-site) %>%
  distinct()

preds$solute <- factor(preds$solute,
                     levels = c("ca",
                                "cl",
                                "mg",
                                "no3",
                                "po4",
                                "so4",
                                "C",
                                "N",
                                "CN"))
levels(preds$solute) <- c("Ca^{`2+`}",
                        "Cl^{`-`}",
                        "Mg^{`2+`}",
                        "NO[3]^{`-`}-N",
                        "PO[4]^{`3-`}-P",
                        "SO[4]^{`2-`}",
                        "'%'*C",
                        "'%'*N",
                        "C:N")


# Tidied model data
tidied <- lmout %>%
  dplyr::select(-data) %>%
  mutate(tid = map(lm4res, tidy)
         ,ano = map(lmres, anova),
         pval = unlist(map(ano, pluck, 4, 2))
         ,gl = map(lm4res, glance)
         ) %>%
  dplyr::select(-lm4res, -preds) %>%
  unnest(gl)

write.csv(tidied, "soils_models/soil_glm_summary_uncorrelated_v2.csv")

# pval text
tidied$pvalt <- ifelse(tidied$pval<0.001,
                       "p<0.001",
                       paste0("p=",round(tidied$pval, 3)))
tidt <- distinct(tidied, solute, pvalt) %>%
  ungroup() %>%
  mutate(x = 0.0,
         y = c(23, 3.1, 8.3, 1.4, 6.1, 5, 35, 2.1, 19))
tidt$solute <- factor(tidt$solute,
                       levels = c("ca",
                                  "cl",
                                  "mg",
                                  "no3",
                                  "po4",
                                  "so4",
                                  "C",
                                  "N",
                                  "CN"))
levels(tidt$solute) <- c("Ca^{`2+`}",
                          "Cl^{`-`}",
                          "Mg^{`2+`}",
                          "NO[3]^{`-`}-N",
                          "PO[4]^{`3-`}-P",
                          "SO[4]^{`2-`}",
                          "'%'*C",
                          "'%'*N",
                          "C:N")

# Plot of predictions
predp <- ggplot(data = preds,
                aes(x = zw,
                    y = conc)) +
  facet_wrap(~solute, scale = "free_y",
             labeller = label_parsed)+
  geom_line() +
  geom_ribbon(aes(ymin=lwr_boot,ymax=upr_boot),
                              alpha=0.05) +
  geom_text(data = tidt,
            aes(x = x,
                y = y,
                label = pvalt),
             size = 3) +
  scale_color_viridis_d() + 
  ylab(expression("Predicted concentration (mg"~L^{-1}*"or -)")) +
  xlab("Relative elevation above water table (m)")
predp

ggsave(predp, 
       filename = "soils_models/soil_predictions_ribbon_v2.tiff",
       device = "tiff",
       dpi = 300)

df$type <- substr(df$site,1,1)

dfp <- df

dfp$solute <- factor(dfp$solute,
                     levels = c("ca",
                                "cl",
                                "mg",
                                "no3",
                                "po4",
                                "so4",
                                "C",
                                "N",
                                "CN"))
levels(dfp$solute) <- c("Ca^{`2+`}",
                         "Cl^{`-`}",
                         "Mg^{`2+`}",
                         "NO[3]^{`-`}-N",
                         "PO[4]^{`3-`}-P",
                         "SO[4]^{`2-`}",
                         "'%'*C",
                         "'%'*N",
                         "C:N")

dfp <- dfp %>%
  group_by(site, solute) %>%
  mutate(conc_s = scale(conc))

conpts2 <- ggplot(dfp, 
                 aes(x=zw,y=conc_s,
                     colour=site,group=site))+
  geom_point() +
  scale_color_viridis_d() + 
  facet_grid(solute~type, scales = "free_y",
             labeller = label_parsed) +
  xlab("Relative elevation (m)") +
  geom_smooth(method = "lm",
              show.legend = FALSE) +
  ylab("Scaled concentration (mg/L or -)")
conpts2
ggsave(conpts2,
       filename = "soils_models/conc_elev_all.tiff",
       device = "tiff",
       height = 7.5,
       width = 6.5,
       units = "in")




dft <- df %>% 
  group_by(site, solute) %>% 
  mutate(dat.med = median(conc))

dft$dat.med.res<-abs(dft$conc-dft$dat.med)

# Then we run an ANOVA, and post-hoc if necessary:
levene.dat.aov<- dft %>%
  mutate(type = substr(site, 1, 1)) %>%
  group_by(solute) %>%
  nest() %>%
  mutate(an = map(data, ~aov(.$dat.med.res~.$type)),
         td = map(an, TukeyHSD),
         t = map(td, tidy))%>%
  unnest(t)

write.csv(levene.dat.aov, "soils_levene.csv")
  aov(dat.med.res~site*solute,dft)
summary(levene.dat.aov)
TukeyHSD(levene.dat.aov)

library(car)

leveneTest(conc~site*solute, data = df)$estimate


