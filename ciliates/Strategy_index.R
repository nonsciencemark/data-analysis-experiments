
# Script to standardize ciliate data by species and make a Strategy index per treatment 
# (strategy index = sum per row of -linearity, -speed, area and aspact ratio)

library(tidyverse)
library(lubridate)
library(mgcv)
library(vegan)

# THE QUESTION IS: CAN PER CAPITA GROWTH RATE BE EXPLAINED BY TRAIT MEANS AND/OR VARIANCES??
# H0: IT CANNOT


# Data entry:
data_cilia <-  read.csv("C:/Users/debruin/OneDrive - UCL/Bio/UCL/GitHub/TeamWork/data-analysis-experiments/ciliates/OverviewTraitsFull.csv") %>%
  filter(ID_spec != "Spiro_5") %>%
  #filter(Treatment == "1") %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(pcgr = log(lead(Parts_permL, 1)/Parts_permL)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Parts_permL) %>% 
  filter(pcgr > -5) %>% #Kick out extremely negative pcgr
  ungroup() %>%
  mutate(Species = as.factor(str_replace_all(str_sub(ID_spec, start = 1, end = 5), "[[:punct:]]", ""))) %>%
  group_by(ID_spec) %>%
  mutate_at(c("mean_speed", "mean_area", "mean_ar", "mean_linearity", "sd_speed", "sd_area", "sd_ar", "sd_linearity", "pcgr", "Parts_permL", "Days_fromstart"), ~ (scale(.) %>% as.vector)) %>%
  ungroup() %>%
  mutate(Strategy_Go = (mean_speed + mean_linearity), 
         Strategy_Stay = (mean_area + mean_ar),
         Strategy_index = (-mean_speed + -mean_linearity + mean_area + mean_ar)) %>% # index by summing speed, linearity, size and shape
  mutate(phase = cut_interval(log10(Parts_permL), n=5, labels = FALSE)) %>%
  select(Species, Atrazine, Temp, Treatment, ID_spec,  Days_fromstart, contains(c("mean_", "sd_")), Strategy_Stay, Strategy_Go, Strategy_index, pcgr, Parts_permL) # without sd for now! All standardized per Species

# convert to factors:
data_cilia$Temp <- as.factor(data_cilia$Temp)
data_cilia$Atrazine <- as.factor(data_cilia$Atrazine)
data_cilia$Treatment <- as.factor(data_cilia$Treatment)
data_cilia$Species <- as.factor(data_cilia$Species)


# Plots:
ggplot(data = data_cilia, mapping = aes(x = Atrazine, y = Strategy_Go, color = ID_spec)) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "pointrange") +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  facet_wrap(facets = .~ Species) +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = data_cilia, mapping = aes(x = Atrazine, y = Strategy_Stay, color = ID_spec)) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "pointrange") +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  facet_wrap(facets = .~ Species) +
  labs() +
  theme_classic(base_size = 24)


ggplot(data = data_cilia, mapping = aes(x = Strategy_Stay, y = Strategy_Go), color = as.factor(Temp)) +
  geom_point(aes(color = Temp) )+  
  facet_wrap(facets = .~ ID_spec) +
  geom_smooth(aes(group = Temp, color = Temp), method = "lm") +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = data_cilia, mapping = aes(x = Atrazine, y = Strategy_index, color = ID_spec)) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "pointrange") +
  facet_wrap(facets = .~ Species + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = subset(data_cilia, Species == "Loxo" | Species == "Para"), mapping = aes(x = pcgr, y = Strategy_Go)) +
  stat_summary(aes(x = pcgr, y = Strategy_Go, shape = Temp, color = Atrazine), 
               fun.data = "mean_se",
               geom = "point", 
               show.legend = TRUE) +
  geom_smooth(aes(group = Atrazine, color = Atrazine), method = "loess", se = F) +
  facet_wrap(facets = .~ ID_spec + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = subset(data_cilia, Species == "Spiro" | Species == "Tetra"), mapping = aes(x = pcgr, y = Strategy_Go)) +
  stat_summary(aes(x = pcgr, y = Strategy_Go, shape = Temp, color = Atrazine), 
               fun.data = "mean_se",
               geom = "point", 
               show.legend = TRUE) +
  geom_smooth(aes(group = Atrazine, color = Atrazine), method = "lm", se = T) +
  facet_wrap(facets = .~ ID_spec + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = data_cilia, mapping = aes(x = pcgr, y = Strategy_Stay)) +
  geom_point(aes(shape = Temp, color = ID_spec), stat = "summary", size = 1, show.legend = TRUE) +
  geom_smooth(aes(group = ID_spec, color = ID_spec), se = F) +
  facet_wrap(facets = .~ Species + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = subset(data_cilia, Species == "Spiro" | Species == "Tetra"), mapping = aes(x = pcgr, y = Strategy_index)) +
  stat_summary(aes(x = pcgr, y = Strategy_index, color = Atrazine), 
               fun.data = "mean_se",
               geom = "point", 
               show.legend = TRUE) +
  geom_smooth(aes(group = Atrazine, color = Atrazine), method = "loess", se = F) +
  facet_wrap(facets = .~ ID_spec + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)

ggplot(data = subset(data_cilia, Species == "Loxo" | Species == "Para"), mapping = aes(x = pcgr, y = Strategy_index)) +
  stat_summary(aes(x = pcgr, y = Strategy_index, color = Atrazine), 
               fun.data = "mean_se",
               geom = "point", 
               show.legend = TRUE) +
  geom_smooth(aes(group = Atrazine, color = Atrazine), method = "loess", se = F) +
  facet_wrap(facets = .~ ID_spec + Temp) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  labs() +
  theme_classic(base_size = 24)
