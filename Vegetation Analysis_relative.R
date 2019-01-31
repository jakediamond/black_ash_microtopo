# 
# Author: Jake Diamond
# Purpose: To analyze veg and chem data continuously with elev
# Date: December 6, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(vegan)
library(lubridate)

# Load veg data
df <- readRDS("veg_data.RDS")
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"

# Load elevation data
elev3 <- readRDS("elevations_wt")

# Get veg data into format for vegan diversity function...no moss
df_div <- df %>%
  filter(moss != 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Get moss data into format for vegan diversity function
df_moss <- df %>%
  filter(moss == 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Function to compute diversity
div_fun <- function(data) {
  badcols <- c("site", "hu.ho", "point")
  mat <- data %>%
    ungroup() %>%
    dplyr::select(-one_of(badcols))
  shann <- diversity(mat)
  simp <- diversity(mat, "simpson")
  rich <- sum(mat > 0)
  x <- data.frame(shannon = shann, 
                  simpson = simp, 
                  richness = rich)
  return(x)
}

# Calculate diversity by site and hummock/hollow
div_site <- df_div %>%
  group_by(site, hu.ho) %>%
  do(div_fun(.))

div_point <- df_div %>%
  group_by(site, point) %>%
  do(div_fun(.)) %>%
  left_join(df %>% 
              dplyr::select(site, point, depth) %>%
              unique()) %>%
  left_join(elev3, by = c("site", "point")) %>%
  ungroup()

div_point$depth2 <- as.numeric(as.character(div_point$depth))
div_point$edge <- ifelse(is.na(div_point$edge), 0, 1)
div_point$log <- ifelse(is.na(div_point$log), 0, 1)
div_point$notree <- ifelse(is.na(div_point$notree), 0, 1)
div_point$intermediate <- ifelse(is.na(div_point$intermediate), 
                                 0, 
                                 1)

# Model
lm_fun <- function(df) {
  lm(richness ~ z_relh_mean, data = df)
}

# Get linear model info
lms <- div_point %>%
  group_by(site) %>%
  nest() %>%
  mutate(model = map(data, lm_fun)) %>%
  mutate(
    adj.r.squared = map_dbl(model, 
                            ~ signif(summary(.)$adj.r.squared, 
                                     3)),
    intercept = map_dbl(model, ~ signif(.$coef[[1]], 3)),
    slope = map_dbl(model, ~ signif(.$coef[[2]], 3)),
    pvalue = map_dbl(model, ~ signif(summary(.)$coef[2, 4], 3)) 
    ) %>%
  dplyr::select(-data, -model)

colnames(lms) <- c("site", "rsq", "b", "m", "pval")
  
lms$ptext <- ifelse(lms$pval < 0.001, 
                          "<0.001", 
                          paste0("=", round(lms$pval, 3)))

# Quick plot
nomoss_elev <- ggplot(div_point, aes(x = z_relh_mean, 
                                     y = richness)) + 
  geom_point(shape = 21, 
             aes(fill = hu.ho,
                 alpha = factor(log))) +
  geom_smooth(method = "lm") + 
  geom_text(data = lms,
    aes(x = 1.4, y = 10, 
        label = paste0("list(y==~", m, "*x", b, ", ",
                       "~italic(R)^2==", rsq,
                       ")")),
    parse = TRUE) +
  geom_text(data = lms,
            aes(x = 1.4, y = 8,
                label = paste0("p ", ptext))) +
  scale_alpha_manual(breaks = c(0, 1),
                    values = c(1, 0.4)) +
  scale_fill_manual(breaks = c("ho", "hu"),
                      values = c("white", "black")) +
  theme_bw() + 
  xlab("Relative elevation (m)") +
  ylab("Richness") +
  theme(legend.position = c(0.75, 0.2)) + 
  facet_wrap(~site, scales = "free_x")

nomoss_elev




# Dataframe for p-values to plot
df_p <- ttests %>%
  select(site, p.value) %>%
  mutate(hu.ho = "hollow",
         y = 1.5)

# Reorder for facets by hydroperiod
levs_hp = c("T3", "D4", "L3", "L1", "D2", "L2",
           "T1", "D3", "T2", "D1")
levs_mean = c("L2", "L1", "T3", "D4", "L3", "D2",
            "T1", "D3", "T2", "D1")
levs_median = c("T3", "L1", "D4", "L3", "D2", "T1",
            "L2", "T2", "D1", "D3")

div_point2$site_f = factor(div_point2$site, 
                         levels = levs_median)

df_p$site_f = factor(df_p$site, 
                           levels = levs_median)

# Plotting Results, no moss for shannon and richness
nomoss <- ggplot(div_point2, aes(x = hu.ho, 
                                 y = shannon, 
                                 fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  geom_text(data = df_p,
            aes(x = hu.ho,
                y = y,
                label = paste0("p = ", round(p.value, 3))),
            show.legend = FALSE) +
  theme_bw() + 
  xlab("") +
  ylab("Shannon Diversity") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = .3) +
  facet_wrap(~site_f, ncol = 5, nrow = 2)

nomoss
ggsave("shannon_sfs.tiff", plot = nomoss,
       device = "tiff",
       dpi = 600)

nomoss2 <- ggplot(div_point, aes(x = hu.ho, y = richness, fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site) + labs(title = "No Moss", 
                           y = "Richness")
ggsave("no_moss_richness.png", plot = nomoss2)

# Plotting Results, no moss and no intermediate points for shannon and richness
nomoss_nointerm <- ggplot(div_point, aes(x = hu.ho, y = shannon, fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site) + labs(title = "No Moss, No Intermediate", 
                           y = "Shannon Diversity")
ggsave("no_moss_no_interm_shannon.png", plot = nomoss_nointerm)

nomoss_nointerm2 <- ggplot(div_point, aes(x = hu.ho, y = richness, fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site) + labs(title = "No Moss, No Intermediate", 
                           y = "Richness")
ggsave("no_moss_no_interm_richness.png", plot = nomoss_nointerm2)

# Plotting Results, moss richness
div_point_moss <- df_moss %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.))

moss <- ggplot(div_point_moss, 
               aes(x = hu.ho, y = shannon, fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site) + labs(title = "Moss", 
                           y = "Shannon Diversity")
ggsave("moss_shannon.png", plot = moss)

moss2 <- ggplot(div_point_moss, aes(x = hu.ho, y = richness, fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site) + labs(title = "Moss", 
                           y = "Richness")
ggsave("moss_richness.png", plot = moss2)



# NMDS Analysis -----------------------------------------------------------
# Get community matrix with sites-hu.ho as index
df_nmds <- df %>%
  dplyr::select(-(10:14)) %>%
  # dplyr::filter(moss != 1) %>%
  left_join(elev3 %>% 
              dplyr::select(site, point, hu.ho, intermediate),
            by = c("site", "point")) %>%
  mutate(hu.ho = ifelse(!is.na(intermediate),
                        "lawn",
                        hu.ho)) %>%
  group_by(site, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0) %>%
  dplyr::rename_all(funs(make.names(.))) %>%
  dplyr::select(-3, -(86:96), -58, -70, -70, -(102:113),-(39:55))
df_nmds<- df_nmds[, colSums(df_nmds != 0) > 0]

groups <- data.frame(site = df_nmds$site,
                     hu_ho = df_nmds$hu.ho)

rownames(df_nmds) <- paste(df_nmds$site, 
                           df_nmds$hu.ho, 
                           sep = ".")
df_nmds$site <- NULL
df_nmds$hu.ho <- NULL

# Run NMDS
nmds <- metaMDS(df_nmds,
                k = 4,
                distance = "bray",
                trymax = 50)

stressplot(nmds)
ordiplot(nmds,type="n")
ordihull(nmds, groups=groups$hu_ho, draw = "polygon", col="grey90", label=F)
orditorp(nmds,display="species",col="red",air=0.01)
orditorp(nmds,display="sites",col=c(rep("green",8),rep("blue",9)),
         air=0.01,cex=1.25)


# Spread of points
plot.new()
ordiellipse(nmds, groups =groups$hu_ho,
                   kind = "se")
dev.off()
#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
#Generate ellipse points
df_ell <- data.frame()
for(g in levels(groups$hu_ho)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, 
                    cbind(as.data.frame(with(groups[groups$hu_ho == g,],
                                             veganCovEllipse(ord[[g]]$cov,
                                                             ord[[g]]$center,
                                                             ord[[g]]$scale))), 
                          hu_ho=g))
  }
}
NMDS.mean=aggregate(nmds[,1:2],list(group=groups$hu_ho),mean)

# > NMDS.mean
# group          x          y
# 1     T -0.2774564 -0.2958445
# 2     V  0.1547353  0.1649902

#Now do the actual plotting
library(ggplot2)

shape_values<-seq(1,11)

p<-ggplot(data=nmds,aes(x,y,colour=hu_ho))
p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
p<-p+geom_point(aes(shape=Depth))+scale_shape_manual(values=shape_values)+theme_bw() 
pdf("NMDS.pdf")
print(p)
dev.off()






# Test for NMDS
# dissimilarities
dis_df <- df_div %>%
  ungroup() %>%
  dplyr::filter(site == "D1") %>%
  dplyr::select(-(1:4))
dis_df <- dis_df[, colSums(dis_df != 0) > 0]
rankindex(scale(dis_df), dis_df)
envfit()

dis <- vegdist(dis_df)
mds0 <- isoMDS(dis)
