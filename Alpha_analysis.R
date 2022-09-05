library(tidyverse)
data <- read.csv("DIVERCE_TdB_Ciliates_Traits.csv") %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("sd_", "ID", "Temp", "Atrazine", "r1", "a11"))) %>%
  pivot_longer(contains(c("sd_")), names_to = "trait", values_to= "sd") %>%
  filter(!grepl('Spiro', ID_spec))

ggplot(data) + 
  theme_classic() + 
  aes(x=log10(sd), y=log10(-a11), 
      col=as_factor(Temp), pch=as_factor(Atrazine)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") #, rows = vars(ID_spec)

