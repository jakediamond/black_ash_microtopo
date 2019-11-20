# 
# Author: Jake Diamond
# Purpose: To analyze vegetation field data categorically
# Date: November 16, 2016
# 

# Set Working Directory
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(vegan)
library(lubridate)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(MASS)
library(lme4)

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
  group_by(site, point, plot, depth, hu.ho, moss) %>%
  summarize(rich = sum(tot > 0, na.rm = T))

# Get elevation data
df_r <- readRDS("elevations_wt") %>%
  dplyr::select(site, point, z_relh_mean) %>%
  right_join(df_r, by = c("site", "point")) %>%
  rename(z = z_relh_mean)
df_r$depth <- as.numeric(df_r$depth)

mean(df_r$rich)
# GLM model
mod <- glmmPQL(rich ~ z, 
               random = ~1|site/plot,
               family = poisson(),
               data = df_r)
mod
summary(mod)
plot(mod)
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(mod,type=c("p","smooth")),
             plot(mod,sqrt(abs(resid(.)))~fitted(.),
                  col=colvec[1],
                  type=c("p","smooth"),
             ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(mod,resid(.,type="pearson")~z,
                  type=c("p","smooth")),
             qqnorm(mod,abline=c(0,1),
                    col=colvec[2]))


mod2 <- glmmPQL(rich ~ moss * z, 
               random = ~1|site,
               family = poisson(),
               data = df_r)
summary(mod2)



mod3 <- glmer(rich ~ z*log(depth) + (1|site),
             family = poisson(),
             data = df_r)
summary(mod3)
plot(mod3)

mod4 <- glmer(rich ~ z*log(depth) + (1|site/plot),
              family = poisson(),
              data = df_r)
summary(mod4)
anova(mod4, mod3)

mod5 <- glmmPQL(rich ~ z,
                random = ~ 1|site,
             family = poisson(),
             data = df_r)


ggplot(data = df_r,
       aes(x = z,
           y = rich,
           color = site))+ geom_point() + 
  stat_smooth(method = "lm",
              se = FALSE) + 
  scale_color_viridis_d() + facet_wrap(~site)

newd <- df_r %>%
  group_by(site) %>%
  mutate(n = n(),
         z = seq(0.5, 2, length.out = n))
pred <- predict(mod5, newdata = newd)
ggplot(data = newd,
       aes(x = z,
           y = rich,
           color = site))+ geom_point() + 
  scale_color_viridis_d()

dfv <- df_r %>%
  group_by(site) %>%
  summarize(varz = sd(z, na.rm = T) / mean(z, na.rm = T),
            varr = sd(rich) / mean(rich))

ggplot(data = dfv,
       aes(x = varz,
           y = varr,
           color = site))+ geom_point() + 
  scale_color_viridis_d()
# Just understory
df_div <- df %>%
  filter(moss != 1) %>%
  group_by(site, point, hu.ho, species) %>%
  summarize(sum = sum(number)) %>%
  spread(key = species, 
         value = sum, 
         fill = 0)
# Same for mosses
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

# Caclulate diversity by site, point, and hummock/hollow
div_both <- df_div_all %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.)) %>%
  left_join(df %>%
              select(site, point, depth)) %>%
  distinct() %>%
  mutate(type = str_sub(site, 1, 1))
# Just vegetation
div <- df_div %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.)) %>%
  left_join(df %>%
              select(site, point, depth)) %>%
  distinct() %>%
  mutate(type = str_sub(site, 1, 1))
div$depth <- as.numeric(as.character(div$depth))

# Same for mosses
div_moss <- df_moss %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.)) %>%
  rename(moss_shan = shannon,
         moss_simp = simpson,
         moss_rich = richness)

# Second dataframe for ease of combining later
div_moss2 <- df_moss %>%
  group_by(site, point, hu.ho) %>%
  do(div_fun(.)) %>%
  mutate(moss = 1,
         type = str_sub(site, 1, 1))

# Combine moss and understory together and save
div_all <- left_join(div, div_moss) %>%
  ungroup()
write_rds(div_all, "diversity_data")

# Long format for plotting moss and no moss together
div_all_l <- div %>%
  mutate(moss = 0) %>%
  full_join(div_moss2) %>%
  select(-depth) %>%
  gather(key = div_type, value = value,
         -site, -hu.ho, -point, -moss)

# T-tests
ttests_both <- div_both %>%
  group_by(site) %>%
  spread(key = hu.ho, 
         value = richness, 
         fill = NA) %>%
  do(tidy(t.test(.$hummock, .$hollow)))
# Just veg
ttests <- div %>%
  group_by(site) %>%
  spread(key = hu.ho, 
         value = shannon, 
         fill = NA) %>%
  do(tidy(t.test(.$hummock, .$hollow)))
# T-tests for moss
ttests_moss <- div_moss2 %>%
  group_by(site) %>%
  spread(key = hu.ho, 
         value = richness, 
         fill = NA) %>%
  do(tidy(t.test(.$hummock, .$hollow)))

# Dissimilarity indices
df_dis_all <- df_div_all %>%
  ungroup() %>%
  select(-point) %>%
  group_by(site, hu.ho) %>%
  summarize_all(funs(sum), na.rm = TRUE)
df_dis_all <- df_dis_all[, colSums(df_dis_all != 0) > 0]

dis_all <- df_dis_all %>%
  ungroup() %>%
  select(-hu.ho) %>%
  group_by(site) %>%
  nest() %>%
  mutate(dis_bray = map(data, vegdist),
         dis_sor = map(data, vegdist, binary = TRUE)) %>%
  select(-data) %>%
  unnest()
# Just veg
df_dis <- df_div %>%
  ungroup() %>%
  select(-point) %>%
  group_by(site, hu.ho) %>%
  summarize_all(funs(sum), na.rm = TRUE)
df_dis <- df_dis[, colSums(df_dis != 0) > 0]

dis <- df_dis %>%
  ungroup() %>%
  select(-hu.ho) %>%
  group_by(site) %>%
  nest() %>%
  mutate(dis_bray = map(data, vegdist),
         dis_sor = map(data, vegdist, binary = TRUE)) %>%
  select(-data) %>%
  unnest()
# Same for moss
df_dis_m <- df_moss %>%
  ungroup() %>%
  select(-point) %>%
  group_by(site, hu.ho) %>%
  summarize_all(funs(sum), na.rm = TRUE)
df_dis_m <- df_dis_m[, colSums(df_dis_m != 0) > 0]

dis_m <- df_dis_m %>%
  ungroup() %>%
  select(-hu.ho) %>%
  group_by(site) %>%
  nest() %>%
  mutate(dis_bray = map(data, vegdist),
         dis_sor = map(data, vegdist, binary = TRUE)) %>%
  select(-data) %>%
  unnest()

# Dataframe for p-values
df_p_all <- ttests_both %>%
  select(site, p.value, estimate1, estimate2) %>%
  mutate(ratio = estimate2 / estimate1,
         hu.ho = "hollow",
         y = 1.5,
         type = str_sub(site, 1, 1)) %>%
  left_join(dis)
# Just veg
df_p <- ttests %>%
  select(site, p.value, estimate1, estimate2) %>%
  mutate(ratio = estimate2 / estimate1,
         hu.ho = "hollow",
         y = 1.5,
         type = str_sub(site, 1, 1)) %>%
  left_join(dis)
# Same for moss
df_p_m <- ttests_moss %>%
  select(site, p.value, estimate1, estimate2) %>%
  mutate(ratio = estimate2 / estimate1,
         hu.ho = "hollow",
         y = 1,
         type = str_sub(site, 1, 1)) %>%
  left_join(dis_m)
# p value text
df_p_all$pvaltext <- ifelse(df_p_all$p.value < 0.001,
                        "p<0.001",
                        paste0("p=", round(df_p_all$p.value, 3)))
df_p$pvaltext <- ifelse(df_p$p.value < 0.001,
                        "p<0.001",
                        paste0("p=", round(df_p$p.value, 3)))

df_p_m$pvaltext <- ifelse(df_p_m$p.value < 0.001,
                        "p<0.001",
                        paste0("p=", round(df_p_m$p.value, 3)))

# Plotting Results, no moss for shannon and richness
p_shan_nm_d <- ggplot(filter(div, type == "D"),
               aes(x = hu.ho, 
                   y = shannon, 
                   fill = hu.ho),
               alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p, type == "D"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p, type == "D"),
            aes(x = hu.ho,
                y = y + 0.1,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory Shannon diversity")

p_shan_nm_l <- ggplot(filter(div, type == "L"),
                      aes(x = hu.ho, 
                          y = shannon, 
                          fill = hu.ho),
                      alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p, type == "L"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p, type == "L"),
            aes(x = hu.ho,
                y = y + 0.1,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory Shannon diversity")

p_shan_nm_t <- ggplot(filter(div, type == "T"),
                      aes(x = hu.ho, 
                          y = shannon, 
                          fill = hu.ho),
                      alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p, type == "T"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p, type == "T"),
            aes(x = hu.ho,
                y = y + 0.1,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory Shannon diversity")

p_shan_nm <- ggdraw() +
  draw_plot(p_shan_nm_d + rremove("x.text") + rremove("x.title") +
              rremove("y.title"), 
            x = 0.036, y = 0.66, width = 0.96, height = 0.33) +
  draw_plot(p_shan_nm_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.355, width = 0.77, height = 0.33) +
  draw_plot(p_shan_nm_t + rremove("y.title"),
            x = 0.038, y = 0, width = 0.7365, height = 0.38)

ggsave(plot = p_shan_nm,
       filename = "shannon_diversity_no_moss.tiff",
       device = "tiff",
       dpi = 300,
       width = 5, 
       height = 4,
       units = "in")

# Plotting Results, no moss for shannon and richness
p_shan_m_d <- ggplot(filter(div_moss2, type == "D"),
                      aes(x = hu.ho, 
                          y = shannon, 
                          fill = hu.ho),
                      alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p_m, type == "D"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_m, type == "D"),
            aes(x = hu.ho,
                y = y,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Moss Shannon diversity")

p_shan_m_l <- ggplot(filter(div_moss2, type == "L"),
                      aes(x = hu.ho, 
                          y = shannon, 
                          fill = hu.ho),
                      alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p_m, type == "L"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_m, type == "L"),
            aes(x = hu.ho,
                y = y,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Moss Shannon diversity")

p_shan_m_t <- ggplot(filter(div_moss2, type == "T"),
                      aes(x = hu.ho, 
                          y = shannon, 
                          fill = hu.ho),
                      alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  geom_text(data = filter(df_p_m, type == "T"),
            aes(x = hu.ho,
                y = y - 0.1,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_m, type == "T"),
            aes(x = hu.ho,
                y = y,
                label = paste0("BC=",
                               round(dis_bray, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Moss Shannon diversity")

p_shan_m <- ggdraw() +
  draw_plot(p_shan_m_d + rremove("x.text") + rremove("x.title") +
              rremove("y.title"), 
            x = 0.036, y = 0.66, width = 0.96, height = 0.33) +
  draw_plot(p_shan_m_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.355, width = 0.77, height = 0.33) +
  draw_plot(p_shan_m_t + rremove("y.title"),
            x = 0.038, y = 0, width = 0.7365, height = 0.38)

ggsave(plot = p_shan_m,
       filename = "shannon_moss.tiff",
       device = "tiff",
       dpi = 300,
       width = 5, 
       height = 4,
       units = "in")

# Both moss and vegetation data in one
p_rich_d <- ggplot(filter(div_both, type == "D"),
                     aes(x = hu.ho, 
                         y = richness, 
                         fill = hu.ho),
                     alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  scale_y_continuous(
                     breaks = seq(0, 12, 3)) + 
  geom_text(data = filter(df_p_all, type == "D"),
            aes(x = hu.ho,
                y = 8,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_all, type == "D"),
            aes(x = hu.ho,
                y = 9,
                label = paste0("SI=",
                               round(dis_sor, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory and moss richness")

p_rich_l <- ggplot(filter(div_both, type == "L"),
                     aes(x = hu.ho, 
                         y = richness, 
                         fill = hu.ho),
                     alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  scale_y_continuous(
                     breaks = seq(0, 12, 3)) + 
  geom_text(data = filter(df_p_all, type == "L"),
            aes(x = hu.ho,
                y = 8,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_all, type == "L"),
            aes(x = hu.ho,
                y = 9,
                label = paste0("SI=",
                               round(dis_sor, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory and moss richness")

p_rich_t <- ggplot(filter(div_both, type == "T"),
                     aes(x = hu.ho, 
                         y = richness, 
                         fill = hu.ho),
                     alpha = 0.8) + 
  geom_bar(position = "dodge", 
           stat = "summary", 
           fun.y = "mean",
           color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = .1) +
  scale_fill_manual(breaks = c("hollow", "hummock"), 
                    values = c("white", "dark grey")) +
  theme_bw() + theme(axis.title.x = element_blank(), 
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_blank(),
                     strip.text.x = element_text(margin =
                                                   margin(0,0,0,0, "cm")),
                     legend.position = "none") +
  scale_y_continuous(
                     breaks = seq(0, 12, 3)) + 
  geom_text(data = filter(df_p_all, type == "T"),
            aes(x = hu.ho,
                y = 8,
                label = pvaltext),
            show.legend = FALSE,
            size = 2) +
  geom_text(data = filter(df_p_all, type == "T"),
            aes(x = hu.ho,
                y = 9,
                label = paste0("SI=",
                               round(dis_sor, 2))),
            show.legend = FALSE,
            size = 2) +
  facet_wrap(~site, ncol = 4) + 
  ylab("Understory and moss richness")

p_rich <- ggdraw() +
  draw_plot(p_rich_d + rremove("x.text") + rremove("x.title") +
              rremove("y.title"), 
            x = 0.036, y = 0.66, width = 0.96, height = 0.33) +
  draw_plot(p_rich_l + rremove("x.text") + rremove("x.title"),
            x = 0, y = 0.355, width = 0.77, height = 0.33) +
  draw_plot(p_rich_t + rremove("y.title"),
            x = 0.038, y = 0, width = 0.7365, height = 0.38)

ggsave(plot = p_rich,
       filename = "richness_all.tiff",
       device = "tiff",
       dpi = 300,
       width = 5, 
       height = 4,
       units = "in")

# NMDS Analysis -----------------------------------------------------------
df_nmds <- df_div_all %>%
  ungroup() %>%
  select(-point) %>%
  group_by(site, hu.ho) %>%
  summarize_all(funs(mean), na.rm = TRUE)

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
                binary = TRUE,
                trymax = 50)

stressplot(nmds)
ordiplot(nmds,type="n")
ordihull(nmds, groups=groups$hu_ho, draw = "polygon", col="grey90", label=F)
orditorp(nmds,display="species",col="red",air=0.01)
orditorp(nmds,display="sites",col=c(rep("green",8),rep("blue",9)),
         air=0.01,cex=1.25)

