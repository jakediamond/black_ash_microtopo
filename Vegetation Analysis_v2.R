# 
# Author: Jake Diamond
# Purpose: To analyze vegetation field data categorically
# Date: November 16, 2016
# 

# Set Working Directory
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(dplyr)
library(tidyr)
library(broom)
library(vegan)
library(ggplot2)
library(lubridate)

# Load veg data
df <- read.csv("veg_data_for_analysis.csv")
df$X <- NULL

# Get data into format for vegan diversity function
df_div <- df %>%
  filter(moss != 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Calculate moss richness
df_moss <- df %>%
  filter(moss == 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(percent)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)

# Function to compute diversity
div_fun <- function(data) {
  # badcols <- c("id", "plot", "position", 
  #              "percent", "depth", "notes", 
  #              "point", "edge", "log", 
  #              "notree", "intermediate", "moss", 
  #              "site", "hu.ho")
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

div_point2 <- df_div %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.)) %>%
  left_join(df %>% 
              select(point, depth)) %>%
  unique()
div_point2$depth2 <- as.numeric(as.character(div_point2$depth))

# T-tests
ttests <- div_point2 %>%
  group_by(site) %>%
  spread(key = hu.ho, 
         value = shannon, 
         fill = NA) %>%
  do(tidy(t.test(.$hummock, .$hollow)))

ggplot(data = div_point, aes(x = depth2, y = richness, color = site)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

# Load water table data to compare relative depths across sites
# wt_df <- read.csv("HydroData/water_table_new_sites.csv")
# wt_df$datetime <- as.POSIXct(wt_df$datetime)
# wt_summary <- wt_df %>%
#   mutate(week = week(datetime)) %>%
#   group_by(site, year, week) %>%
#   summarize(meanwt = mean(waterlevel, na.rm = TRUE)) %>%
#   group_by(site, year) %>%
#   summarize(hydroperiod = sum(meanwt > 0, na.rm = TRUE),
#             sd = sd(meanwt, na.rm = TRUE))

# Order sites by hydroperiod, low-to-high
div_point2$site = factor(div_point2$site, 
                         levels = c("T3", "B1", "D4", "B6", "D2", "B3",
                                    "T1", "D1", "D3", "T2")
                         )


# Plotting Results, no moss for shannon and richness
nomoss <- ggplot(div_point2, aes(x = hu.ho, 
                                 y = shannon, 
                                 fill = hu.ho)) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.title.y = element_text(face = "bold", 
                                                 vjust = 0.6),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(size = 1,
                                                       linetype = "dashed"),
                     strip.text = element_text(face = "bold"),
                     legend.position = "none") +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("black", "grey")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .1) +
  facet_wrap(~site, ncol = 5, nrow = 2) + labs(y = "Shannon Diversity")
ggsave("shannon_poster.png", plot = nomoss)

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

moss <- ggplot(div_point_moss, aes(x = hu.ho, y = shannon, fill = hu.ho)) + 
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
