library(tidyverse)
library(lubridate)
library(mgcv)
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# CILIATES ------
## Import data and compute pcgr----
data_cilia <- read.csv("ciliates/OverviewTraitsFull.csv") %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(pcgr = log(lead(Count, 1)/Count)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Count) %>% 
  ungroup() %>%
  group_by(ID_spec) %>%
  mutate(phase = cut_interval(log10(Count), n=5, labels = FALSE)) %>%
  filter(pcgr > -5) %>% #Kick out extremely negative pcgr
  pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  separate(trait, into=c("stat", "trait"), sep="_") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean)

ggplot(data_cilia) +
  scale_shape_manual(values = 0:8) +
  theme_bw() + 
  geom_point(aes(x=mean, y=pcgr, pch=as_factor(Treatment),
                 col=as_factor(phase)), alpha=0.5) + 
  geom_smooth(method=lm, aes(x=mean, y=pcgr, 
                             colour=as_factor(phase)), 
              formula=y ~ poly(x,1)) + #, se=F
  facet_wrap(vars(trait), scales="free") #, rows = vars(ID_spec)strain ~ 

ggplot(data_cilia) +
  scale_shape_manual(values = 0:8) +
  theme_bw() + 
  geom_point(aes(x=cv, y=pcgr, pch=as_factor(Treatment),
                 col=as_factor(phase)), alpha=0.5) + 
  geom_smooth(method=lm, aes(x=cv, y=pcgr, colour=as_factor(phase)), 
              formula=y ~ poly(x,1), se=F) +
  facet_wrap(vars(trait), scales="free")

# CYANO ------
## Import data and compute pcgr----
data_cyano <- read.csv("cyanobacteria/mono_data.csv") %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  group_by(strain, treat) %>%
  mutate(day = day - min(day)) %>%
  mutate(pcgr = log(lead(population.mean,1)/population.mean)/(lead(day,1)-day),
         .after = population.mean) %>%
  ungroup() %>%
  group_by(strain) %>%
  mutate(phase = cut_interval(log10(population.mean), n=5, labels = FALSE)) %>%
  select(-population.sd) %>%
  pivot_longer(FSC.HLin.mean:NIR.R.HLin.sd, names_to="trait") %>%
  separate(trait, into=c("trait", "stat"), sep="in.") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean)

ggplot(data_cyano) +
  theme_bw() + 
  geom_point(aes(x=mean, y=pcgr, pch=as_factor(treat),
                 col=as_factor(phase)), alpha=0.5) + 
  geom_smooth(method=lm, aes(x=mean, y=pcgr, 
                             colour=as_factor(phase)), 
              formula=y ~ poly(x,1)) + #, se=F
  facet_wrap(vars(trait), scales="free") #, rows = vars(ID_spec)strain ~ 

ggplot(data_cyano) +
  theme_bw() + 
  geom_point(aes(x=cv, y=pcgr, pch=as_factor(treat),
                 col=as_factor(phase)), alpha=0.5) + 
  geom_smooth(method=lm, aes(x=cv, y=pcgr, 
                             colour=as_factor(phase)), 
              formula=y ~ poly(x,1), se=F) +
  facet_wrap(vars(trait), scales="free") 
