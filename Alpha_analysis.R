library(tidyverse)
library(lubridate)
data_cilia <- read.csv("ciliates/DIVERCE_TdB_Ciliates_Traits.csv") %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("sd_", "ID", "Temp", "Atrazine", "r1", "a11"))) %>%
  pivot_longer(contains(c("sd_")), names_to = "trait", values_to= "sd") %>%
  filter(!grepl('Spiro', ID_spec))

ggplot(data_cilia) + 
  theme_classic() + 
  aes(x=log10(sd), y=log10(-a11), 
      col=as_factor(Temp), pch=as_factor(Atrazine)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") #, rows = vars(ID_spec)

data_cyano <- read.csv("cyanobacteria/mono_data.csv") %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  group_by(strain, treat) %>%
  mutate(day = day - min(day)) %>%
  
  