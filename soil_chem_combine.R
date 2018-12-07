# 
# Author: Jake Diamond
# Purpose: Combine all soil chem data (CN + extractions)
# Date: December 6, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)

# Load data
df_cn <- read.csv("Soils/CN.csv")
hu <- read.csv("hu.ho2.csv")
df <- read.csv("Soils/meta_for_r.csv", stringsAsFactors = FALSE)
df2 <- read.csv("soil_chem_r.csv")

# Change point column in df to match elev
df$point <- paste(df$plot, df$point, sep = ".")

# Change B sites to L sites
df[df$site == "B1", "site"] <- "L1"
df[df$site == "B3", "site"] <- "L2"
df[df$site == "B6", "site"] <- "L3"

# Rename points to match
soil$Name <- paste(substr(soil$Name, 1, 2), 
                   ".", 
                   substring(soil$Name, 4),
                   sep = "")
# Clean data
soil <- subset(soil, Weight1 > 200)
soil <- subset(soil, !(Name %in% c("je.f soil", 
                                   "lo. std",
                                   "bl.nk")))
soil$Name <- sapply(strsplit(soil$Name, " ", fixed = TRUE),
                    function(x) (x[1]))
colnames(soil)[2] <- "point"

# Average reps
soil_avg <- soil %>%
  group_by(point) %>%
  summarize(n = mean(X..N, na.rm = TRUE),
            c = mean(X..C, na.rm = TRUE),
            cn = mean(CNRatio, na.rm = TRUE))

# Combine data
df <- left_join(soil_avg, hu)
df <- subset(df, !is.na(site))

# Plot data
cn <- ggplot(df %>%
               dplyr::filter(is.na(log)),
             aes(x = hu.ho, 
                 y = c, 
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
  facet_wrap(~site, ncol = 5, nrow = 2) + labs(y = "CN")
cn

ggsave("cn.png", plot = cn)


# Global plot
cn_g <- ggplot(df %>%
                 dplyr::filter(is.na(log)),
               aes(x = hu.ho, 
                   y = n, 
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
  labs(y = "CN")
cn_g
ggsave("cn_global.png", plot = cn_g)


# Analyze data
ttests <- df %>%
  dplyr::filter(is.na(log)) %>%
  spread(key = hu.ho, 
         value = n, 
         fill = NA) %>%
  do(tidy(t.test(.$hummock, .$hollow)))


# Find replicates
df2$sample_id <- as.character(df2$sample_id)
df2$reps <- unlist(strsplit(df2$sample_id, "R"))

# Make calcium numeric because it thinks its a character
df2$ca <- as.numeric(df2$ca)

# Find average and rpd of replicates
df_rep <- df2 %>%
  group_by(reps) %>%
  filter(n() > 1) %>%
  group_by(reps) %>%
  summarize_at(vars(cl:po4),
               funs(mean, diff)) %>%
  gather(key, value, -reps) %>%
  extract(key, c("solute", "calc"), "(^[^_]+)_(.*)$") %>%
  spread(calc, value) %>%
  mutate(rpd = abs(diff) * 100 / mean) %>%
  mutate(reps = as.numeric(reps))

df_rep2 <- df %>%
  rename(reps = sample_id) %>%
  right_join(df_rep)

# Graph of rpd
p_rep <- ggplot(data = df_rep2, 
                aes(x = solute,
                    y = rpd,
                    fill = hu.ho)) +
  geom_bar(stat = "summary", 
           fun.y = "mean",
           position = "dodge") + 
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               width = 0.4,
               position = "dodge") +
  theme_bw() + 
  facet_wrap(~site) +
  xlab("Solute") + 
  ylab("Relative percent difference (%)")

p_rep
# Save plot
ggsave(plot = p_rep, 
       "rpd_for_replicates_huho_site.tiff",
       device = "tiff",
       width = 8,
       height = 6,
       units = "in")

# Get all data in one place
df3 <- df2 %>%
  mutate(reps = as.numeric(reps)) %>%
  group_by(reps) %>%
  filter(n() == 1) %>%
  gather(solute, mean, -reps, -sample_id) %>%
  full_join(df_rep) %>%
  select(-rpd, -diff) %>%
  mutate(sample_id = reps,
         sample_id = as.integer(sample_id))

df <- df3 %>%
  spread(solute, mean) %>%
  right_join(df)

df_l <- df %>%
  gather(key = solute, 
         value = conc,
         ca, cl, mg, no3, po4, so4)

# Strip panel names for plot
levels(df_l$solute) <- c("Ca^{`2+`}",
                         "Cl^{`-`}",
                         "Mg^{`2+`}",
                         "NO[3]^{`-`}",
                         "PO[4]^{`3-`}",
                         "SO[4]^{`2-`}")