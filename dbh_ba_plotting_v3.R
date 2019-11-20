# 
# Purpose: Analyze and plot dbh as a function of elevation/hummock or hollow
# Author: Jake Diamond
# Date: February 27, 2019
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
# setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")
setwd("C:/Users/jake.diamond/Dropbox/Projects/EAB/Data")

# Load Libraries
library(Hmisc)
library(broom)
library(viridis)
library(ggrepel)
library(scales)
library(tidyverse)
# Set plot theme
theme_set(theme_bw()+
            theme(panel.grid = element_blank(),
                  strip.text.x = element_text(margin = margin(0,0,0,0, "cm"))))

# Load xyz/hummock data of dbh by site
df <- read_rds("dbh_hummock_match_minz")
humsites <- unique(df$site)

# Check on hummock hollow data...use cleaned data
huclean <- read_csv("hummock_stats_clean.csv") %>%
  filter(site %in% c("B1", "B3")) %>%
  mutate(site = ifelse(site == "B1",
                       "L1",
                       ifelse(site == "B3",
                              "L2",
                              site)))

# Load site level basal area (m2/ha) data
ba <- read_csv("Forestry/forest_hydro_summary.csv")

# Get all elevation relative to mean wt
df <- df %>%
  left_join(ba %>%
              dplyr::select(site, mean)) %>%
  mutate(zrel = zuse - mean)

# Quick ba summary
bas <- ba %>%
  rename(tph = "trees_per_ha") %>%
  group_by(site) %>%
  dplyr::select(overstory2, midstory2, median, mean, tph) %>%
  rename(Canopy = overstory2, Midstory = midstory2) %>%
  gather(story, ba, -median, -site, -mean, -tph)
bas$type <- str_sub(bas$site, 1, 1)

# Linear models of median vs basal area
lin_mods <- bas %>%
  group_by(story) %>%
  nest() %>%
  mutate(model = map(data, ~lm(ba ~ median, data = .)),
         glance_lm = map(model, glance),  
         rsq = map_dbl(glance_lm, "r.squared"),
         pval = map_dbl(glance_lm, "p.value"),
         tidied = map(model, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE) %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate)

# Pvalue text for plot
lin_mods$pval_text <- ifelse(lin_mods$pval < 0.001, 
                            "p<0.001", 
                            paste0("p=", round(lin_mods$pval, 3)))
lin_mods$x <- 0.075
lin_mods$y <- c(50, 5.5)
# Anova among types
mods <- bas %>%
  group_by(story) %>%
  nest() %>%
  mutate(m = map(data, ~lm(ba ~ type, data = .)),
         a = map(m, aov),
         thsd = map(a, TukeyHSD),
         t = map(thsd, tidy)) %>%
  unnest(t)
# Text for graph
ps <- data.frame(type = c("D", "L", "T"),
                 y = 50,
                 label = c("a", "b", "c"))
hs <- data.frame(type = c("D", "L", "T"),
                 y = 10,
                 label = c(0.01, -0.32, -0.04))

typeba <- ggplot(data = dplyr::filter(bas, story == "Canopy"),
       aes(x = type,
           y = ba,
           fill = type)) + 
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",
               color = "black",
               alpha = 0.7) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",
               width = 0.1) +
  geom_text(data = ps,
            aes(x = type,
                y = y,
                label = label),
            color = "black") + 
  geom_text(data = hs,
            aes(x = type,
                y = y,
                label = paste0("bar(WT)==", label, "~m")),
            color = "black",
            parse = TRUE) + 
  scale_fill_viridis(discrete = TRUE) +
  ylab(expression("Canopy basal area ("*m^2*ha^{-1}*")")) +
  xlab("") + 
  theme(legend.position = "none")
typeba
ggsave(typeba,
       filename = "Forestry/type_ba.tiff",
       device = "tiff",
       dpi = 300)

modstph <- bas %>%
  group_by(story) %>%
  nest() %>%
  mutate(m = map(data, ~lm(tph ~ type, data = .)),
         a = map(m, aov),
         thsd = map(a, TukeyHSD),
         t = map(thsd, tidy)) %>%
  unnest(t)

# Linear model of organic layer 
ba$org_depth <- ifelse(is.na(ba$org_depth), 0 ,ba$org_depth)
ba$org_depth <- ba$org_depth / 100
summary(lm(overstory2 ~ org_depth, data = ba))

# Plot of site level basal area vs. water table
bawt <- ggplot(data = bas,
       aes(x = median,
           y = ba)) +
  geom_point(
    aes(color = type),
             size = 2) +
  geom_text_repel(aes(label = site),
                   nudge_x = 0.01) + 
  facet_wrap(~story, scales = "free_y") +
  geom_smooth(data = dplyr::filter(bas, site != "T1"),
              method = "lm", alpha = 0.1) + 
  scale_color_viridis(discrete = TRUE,
                      name = "Site type") +
  theme(legend.position = c(0.08, 0.20),
        legend.background = element_rect(fill = "grey90",
                                         color = "black"),
        legend.margin = margin(0.08,0.1,0.08,0.1, "cm")) +
    geom_text(data = lin_mods,
            aes(x = x,
                y = y,
                label = pval_text),
            show.legend = FALSE,
            size = 4) +
  ylab(expression("Basal area ("*m^2*ha^{-1}*")")) + 
  xlab("Median daily water table (m)")
bawt
ggsave(bawt,
       filename = "Forestry/ba_wt_new_v2.tiff",
       device = "tiff",
       dpi = 300,
       width = 5.5,
       height = 3.5,
       units = "in")

bas_wt <- bas %>%
  group_by(type) %>%
  dplyr::summarize(med = mean(median),
                   medsd = sd(median),
                   mea = mean(mean),
                   measd = sd(mean))
bas2 <- bas_wt %>%
  right_join(bas %>%
               dplyr::select(type, ba, story))

# Group by site type
bawts <- ggplot(data = bas2,
                aes(x = med,
                    y = ba,
                    group = type,
                    color = type)) +
  # geom_boxplot() +
  stat_summary(fun.data = mean_cl_boot) +
  # geom_errorbarh(data = bas_wt,
  #                aes(x = med,
  #                    xmin = med - medsd,
  #                    xmax = med + medsd)) +
  # facet_wrap(~story, scales = "free_y") +
  # geom_smooth(method = "lm", alpha = 0.1) + 
  # scale_y_continuous(trans = log2_trans(),
  #                    breaks = trans_breaks("log2", function(x) 2^x),
  #                    labels = trans_format("log2", math_format(2^.x))) +
  scale_color_viridis(discrete = TRUE,
                      name = "Site type") +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey90",
                                         color = "black"),
        legend.margin = margin(0.08,0.1,0.08,0.1, "cm")) +
  ylab(expression("Basal area ("*m^2*ha^{-1}*")")) + 
  xlab("Median daily water table (m)")
bawts

# Summarize tls dbh by hummock or hollow
df <- df %>%
  left_join(huclean %>%
              dplyr::select(site, id) %>%
              mutate(clean = 1)) %>%
  mutate(
         clean = ifelse(!(site %in% c("L1", "L2")),
                            1,
                            clean)
         , huho = ifelse(is.na(zmax_raw) | is.na(clean),
                        "hollow",
                        "hummock")
         )
# ttests
tts <- df %>%
  dplyr::select(site, huho, dbh) %>%
  group_by(site) %>%
  mutate(i = row_number()) %>%
  spread(huho, dbh) %>%
  dplyr::filter(site != "L2") %>%
  nest() %>%
  mutate(tt = map(data, ~t.test(.$hollow, .$hummock, na.rm = TRUE)),
         t = map(tt, tidy)) %>%
  unnest(t) %>%
  mutate(huho = "hollow")

# Summary histogram of DBH in hollows vs hummocks
dbhuho <- ggplot(data = df,
       aes(x = huho,
           y = dbh,
           fill = huho)) + 
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",
               color = "black") +
  stat_summary(geom = "errorbar", 
               fun.data = mean_se, 
               position = "dodge", width = 0.1) +
  scale_fill_manual(breaks = c("hollow", "hummock"),
                    values = c("white", "dark grey")) +
  facet_wrap(~site) + 
  geom_text(data = tts,
            aes(x = huho,
                y = 30,
                label = paste0("p=", round(p.value, 3)))) +
  ylab("Tree diameter at breast height (cm)") + 
  xlab("") + 
  theme(legend.position = "none")
dbhuho
ggsave(dbhuho,
       filename = "Forestry/dbh_huho_cat.tiff",
       device = "tiff",
       dpi = 300,
       width = 5, 
       height = 4,
       units = "in")

# summary of percent of trees in hollows
dfsum <- df %>%
  group_by(site) %>%
  dplyr::summarize(n = n(),
            hu = sum(huho == "hummock"),
            hurat = hu / n)

# Plot of water table vs ratio of hummock vs hollow
d <- left_join(ba, dfsum)
ggplot(data = d,
       aes(x = median,
           y = hurat)) + geom_text(aes(label = site))

# Compare rates of big and small trees on hu ho
dfsum <- df %>% 
  mutate(clas = ifelse(dbh <= 20,
                "less",
                "more")) %>%
  group_by(site, clas) %>%
  dplyr::summarize(hurat = sum(huho == "hummock") / n()) %>%
  bind_rows(dfsum %>%
              mutate(clas = "all") %>%
              dplyr::select(-n, -hu))

# Histogram of dbh by hummock/hollow and site
histhuho <- ggplot(data = df) + 
  geom_histogram(aes(x = dbh,
                     fill = huho)) + 
  facet_wrap(~site) + 
  geom_text(data = dplyr::filter(dfsum, clas == "all"),
            aes(x = 40,
                y = 20,
                label = paste0("f[hum]==",
                               round(hurat, 2))),
            parse = TRUE,
            size = 3) + 
  geom_text(data = dplyr::filter(dfsum, clas == "less"),
            aes(x = 40,
                y = 18,
                label = paste0("f[hum*`,<20`]==",
                               round(hurat, 2))),
            parse = TRUE,
            size = 3) + 
  geom_text(data = dplyr::filter(dfsum, clas == "more"),
            aes(x = 40,
                y = 16,
                label = paste0("f[hum*`,>20`]==",
                               round(hurat, 2))),
            parse = TRUE,
            size = 3) + 
  scale_fill_manual(name = "Microsite",
                     values = c("grey50", "black"),
                     labels = c("Hollow",
                                "Hummock")) + 
  theme(legend.position = c(0.09, 0.89),
        legend.background = element_rect(fill = "grey90",
                                         color = "black"),
        legend.margin = margin(0.08,0.1,0.08,0.1, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.7, "line")) +
  ylab("Count") + xlab("Tree diameter at breast height (cm)")
histhuho
ggsave(histhuho,
       filename = "Forestry/dbh_huho_bw.tiff",
       device = "tiff",
       dpi = 300,
       width = 5.5,
       height = 4, 
       units = "in")

library(boot)
df_sum <- df %>%
  group_by(site) %>%
  dplyr::summarise(site_area_all = mean(site_area, na.rm = T)) %>%
  ungroup() %>%
  right_join(df) %>%
  select(site, huho, dbh, site_area_all) %>%
  mutate(area = dbh * dbh * 3.14 / 4)

x <- df_sum %>%
  group_by(site, huho, site_area_all) %>%
  nest() %>%
  mutate(b = map(data, ~boot(.$area, sum, R = 500)),
         bci = map(b, boot.ci, type = "norm"),
         statistic = purrr::map(.x = bci,
                                ~ .x$t0), 
         lower_ci = purrr::map(.x = bci,
                               ~ .x$normal[[2]]), 
         upper_ci = purrr::map(.x = bci,
                               ~ .x$normal[[3]])) %>% 
  select(-data, -b, -bci) %>%
  unnest() %>%
  mutate(stat_norm = statistic / site_area_all,
         lower_norm = lower_ci / site_area_all,
         upper_norm = upper_ci / site_area_all,
         diff = statistic - lower_ci,
         diffnorm = stat_norm - lower_norm)


ggplot(data = x,
                 aes(x = huho,
                     y = stat_norm,
                     fill = huho)) + 
  stat_summary(geom = "bar", fun.y = sum, position = "dodge",
               color = "black") +
  geom_errorbar(aes(ymin = lower_norm,
                    ymax = upper_norm),
               position = "dodge", width = 0.1) +
  scale_fill_manual(breaks = c("hollow", "hummock"),
                    values = c("white", "dark grey")) +
  facet_wrap(~site) + 
  # geom_text(data = tts,
  #           aes(x = huho,
  #               y = 30,
  #               label = paste0("p=", round(p.value, 3)))) +
  ylab("Cumulative basal area (m^2 / m^2)") + 
  xlab("") + 
  theme(legend.position = "none")





# look at hum vs hol
huhocomp <- df %>%
  group_by(site) %>%
  dplyr::summarise(site_area_all = mean(site_area, na.rm = T)) %>%
  ungroup() %>%
  right_join(df) %>%
  group_by(site, huho) %>%
  dplyr::summarize(mean = mean(dbh, na.rm = T),
            n = n(),
            se = sd(dbh, na.rm = T) / sqrt(n),
            sum = sum(dbh * dbh * 3.14 / 4, na.rm = T),
            site_area = mean(site_area_all),
            sum_norm = sum / site_area)

# t-tests for hummock-hollow dbh dif (L2 has no hummocks)
ttests_dbh <- df %>%
  group_by(site) %>%
  filter(site != "L2") %>%
  nest() %>%
  mutate(data = map(data, spread, huho, dbh)) %>%
  mutate(
    ttest = map(data, ~ t.test(.x$hummock, .x$hollow)),
    tidied = map(ttest, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE)

# Dataframe for p-values and estimates
df_text <- ttests_dbh %>%
  dplyr::select(site, p.value, estimate)
df_text$p.value <- ifelse(df_text$p.value < 0.001, 
                          "p<0.001", 
                          paste0("p=",round(df_text$p.value, 3)))
df_text$y <- 50
df_text$x <- "hollow"

df2 <- df %>%
  mutate(huho = ifelse(is.na(zmax_raw),
                       "hollow",
                       "hummock"))
# Quick plot of dbh differences
p_dbh <- ggplot(data = df2, 
                aes(x = huho,
                    y = dbh)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~site) +
  geom_text(data = df_text,
            aes(x = x,
                y = y,
                label = paste0("p = ", p.value)),
            show.legend = FALSE) +
  geom_text(data = huhocomp,
            aes(x = huho,
                y = 44,
                label = paste("list(Delta*z ==",
                              round(estimate, digits=2), ")")),
            show.legend = FALSE,
            parse = TRUE) +
  ylab("DBH (cm)") +
  xlab("")
p_dbh
ggsave(plot = p_dbh,
       filename = "hummock_hollow_dbh.tiff",
       dpi = 300,
       device = "tiff")


# tests
df <- df %>%
  mutate(
    zuse = ifelse(is.na(zmean_raw),
                  z,
                  zmean_raw),
    zuse = zuse - zwell,
    zrel = zuse - mean) 
# Quick edit
dfp <- df %>%
  mutate(huho = ifelse(zuse > 0.4 & !(site %in% c("L1", "L2")),
                       "hummock",
                       huho))

hh <- df %>%
  group_by(site) %>%
  dplyr::filter(huho == "hummock") %>%
  dplyr::summarize(m = median(zrel))

# Quick plot
dbhp <- ggplot(dfp, aes(x = zrel, y = dbh)) + 
  scale_shape_manual(name = "Microsite",
                     breaks = c("hollow", "hummock"),
                     labels = c("Hollow", "Hummock"),
                     values = c(1, 16)) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.background = element_rect(color = "black",
                                         fill = "grey90")) +
  geom_point(aes(shape = huho)) + 
  facet_wrap(~site, scales = "free_x") +
  geom_smooth(method = "lm",
              show.legend = FALSE) +
  ylab("DBH (cm)") + 
  xlab("Relative elevation above water table (m)")
dbhp
ggsave(dbhp,
       filename = "Forestry/dbh_vs_elev_huho.tiff",
       device = "tiff",
       dpi = 300,
       width = 4,
       height = 3,
       units = "in")

mods_z <- df %>%
  group_by(site) %>%
  nest() %>%
  mutate(l = map(data, ~lm(dbh ~ zmaxn, data = .)),
         t = map(l, tidy)) %>%
  unnest(t)



