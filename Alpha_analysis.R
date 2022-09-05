library(tidyverse)
library(lubridate)
# CILIATES ------
data_cilia <- read.csv("ciliates/DIVERCE_TdB_Ciliates_Traits.csv") %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("sd_", "ID", "Temp", "Atrazine", "r1", "a11"))) %>%
  pivot_longer(contains(c("sd_")), names_to = "trait", values_to= "sd") %>%
  filter(!grepl('Spiro', ID_spec))

ggplot(data_cilia) + 
  theme_classic() + 
  aes(x=log10(sd), y=a11, 
      col=as_factor(Temp), pch=as_factor(Atrazine)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") + #, rows = vars(ID_spec)
  geom_smooth(method="lm")

# CYANO ------
data_cyano <- read.csv("cyanobacteria/mono_data.csv") %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  group_by(strain, treat) %>%
  mutate(day = day - min(day)) %>%
  mutate(day.0 = lag(day, 1), .after = day) %>%
  mutate(population.mean.0 = lag(population.mean, 1), 
         .after = population.mean) %>%
  mutate(pcgr = log(population.mean/population.mean.0)/(day-day.0),
         .after = population.mean.0)

#kick out: 
#all pcgr that are too negative (crashing o/t culture) and 
#all data points where log10(population.mean) is consistently low (<3.5) because of lags
data_cyano <- data_cyano %>%
  filter(pcgr>-0.1) %>%
  mutate(lag.phase = ((log10(population.mean.0)<3.5)&(log10(population.mean)<3.5))) %>%
  filter(lag.phase == FALSE)

ggplot(data_cyano) +
  theme_classic() + 
  aes(x=day, y=log10(population.mean), 
      col=as_factor(treat), pch=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(strain), scales="free") #, rows = vars(ID_spec)

ggplot(data_cyano) +
  theme_classic() + 
  aes(x=population.mean.0, y=pcgr, 
      col=as_factor(treat), pch=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(strain), scales="free") #, rows = vars(ID_spec)

#Now compute the alpha's
data_cyano <- data_cyano %>%
  group_by(strain, treat) %>%
  mutate(alpha = -cov(population.mean.0, pcgr) / var(population.mean.0)) %>% #alphas
  group_by(strain, treat) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("strain", "treat", "alpha", ".sd"))) %>%
  select(-contains("population")) %>%
  pivot_longer(contains(c(".sd")), names_to = "trait", values_to= "sd")

ggplot(data_cyano %>% filter(treat %in% c("C", "T"))) + 
  theme_classic() + 
  aes(x=log10(sd), y=log10(alpha), 
      col=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") + 
  geom_smooth(method="lm") #, rows = vars(ID_spec)



