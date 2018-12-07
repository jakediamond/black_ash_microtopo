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
    select(-one_of(badcols))
  shann <- diversity(mat)
  simp <- diversity(mat, "simpson")
  rich <- sum(mat > 0)
  x <- data.frame(shannon = shann, simpson = simp, richness = rich)
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
  left_join(elev, by = c("site", "point")) %>%
  ungroup()

div_point$depth2 <- as.numeric(as.character(div_point$depth))
div_point$edge <- ifelse(is.na(div_point$edge), 0, 1)
div_point$log <- ifelse(is.na(div_point$log), 0, 1)
div_point$notree <- ifelse(is.na(div_point$notree), 0, 1)
div_point$intermediate <- ifelse(is.na(div_point$intermediate), 0, 1)

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
    adj.r.squared = map_dbl(model, ~ signif(summary(.)$adj.r.squared, 3)),
    intercept = map_dbl(model, ~ signif(.$coef[[1]], 3)),
    slope = map_dbl(model, ~ signif(.$coef[[2]], 3)),
    pvalue = map_dbl(model, ~ signif(summary(.)$coef[2, 4], 3)) 
    ) %>%
  select(-data, -model)

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




# Get community matrix with sites-hu.ho as index
df_nmds <- df %>%
  filter(moss != 1) %>%
  group_by(site, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

rownames(df_nmds) <- paste(df_nmds$hu.ho, df_nmds$site, sep = ".")
df_nmds$site <- NULL
df_nmds$hu.ho <- NULL
df_nmds2 <- df_nmds[order(row.names(df_nmds)), ]

example_NMDS <- metaMDS(df_nmds2, # Our community-by-species matrix
                     k=3)
stressplot(example_NMDS)
plot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)
huho <- c(rep("Hollow",10),rep("Hummock",10))
ordiplot(example_NMDS,type="n")
ordihull(example_NMDS, groups=huho, draw = "polygon", col="grey90", label=F)
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
