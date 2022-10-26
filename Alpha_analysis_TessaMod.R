library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# CILIATES ------
## Import data and compute pcgr, and only keep specific species-trait combinations----
data_cilia <- data_cilia %>%
  rename(strain = ID_spec) %>%
  group_by(Atrazine, Temp, strain) %>%
  mutate(pcgr = log(lead(Parts_permL, 1)/Parts_permL)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Parts_permL) %>% 
  filter(pcgr > -5) %>% #Kick out extremely negative pcgr
  ungroup() %>%
  separate(strain, into = c("species", "sp.strain"), remove=F) %>%#link to species
  group_by(strain) %>%
  mutate(phase = cut_interval(log10(Parts_permL), n=5, labels = FALSE)) %>%
  pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  separate(trait, into=c("stat", "trait"), sep="_") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean)%>%
  filter(((species=="Spiro")&(trait=="linearity"))|((species=="Tetra")&(trait=="area"))
         |((species=="Para")&(trait=="speed"))|((species=="Loxo")&(trait=="linearity")))

## fit a reference model: of dd for C data only, and add to the original data frame:------
ref_model <- lm(data_cilia %>% filter(Treatment==1), 
                formula = pcgr ~ poly(Parts_permL,1)*strain + strain)
data_cilia$pcgr_ref <- predict.lm(ref_model, newdata = data_cilia)
data_cilia <- data_cilia %>%
  mutate(delta = pcgr - pcgr_ref)
#plot the data and overlay the reference model of dd
ggplot(data_cilia) +
  theme_bw() + 
  geom_point(aes(x=Parts_permL, y=pcgr)) + 
  geom_smooth(method=lm, aes(x=Parts_permL, y=pcgr), 
              formula=y ~ poly(x,1), se=F) + #
  geom_line(aes(x=Parts_permL, y=pcgr_ref)) +
  facet_wrap(vars(Treatment, strain), scales="free")

#now check if the deviation from the ref model (delta) 
#depends on the trait value
ggplot(data_cilia, mapping = aes(x=Strategy_index, y=delta, col=strain)) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  geom_point(aes(pch=as_factor(Temp)), stat = "summary") + 
  geom_smooth(method=lm, aes(x=mean, y=delta, col=as_factor(Treatment), group = strain),
              formula=y ~ poly(x,1), se=F) + #
  facet_wrap(vars(Atrazine)) #, rows = vars(ID_spec)strain ~ 
#It doesn't 

# CYANO ------
## Import data and compute pcgr, and only keep specific species-trait combinations----
data_cyano <- read.csv("cyanobacteria/mono_data.csv") %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  group_by(strain, treat) %>%
  mutate(day = day - min(day)) %>%
  mutate(pcgr = log(lead(population.mean,1)/population.mean)/(lead(day,1)-day),
         .after = population.mean) %>%
  filter(!is.na(pcgr)) %>%
  ungroup() %>%
  mutate(species = case_when(strain%in%c(2375,2524)~"one",
                             TRUE ~ "two")) %>%#link to species
  mutate(strain = as_factor(strain)) %>%
  group_by(strain) %>%
  mutate(phase = cut_interval(log10(population.mean), n=5, labels = FALSE)) %>%
  select(-population.sd) %>%
  pivot_longer(FSC.HLin.mean:NIR.R.HLin.sd, names_to="trait") %>%
  separate(trait, into=c("trait", "stat"), sep="in.") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean) %>%
  filter(((species=="one")&(trait=="RED.R.HL"))|((species=="two")&(trait=="YEL.B.HL")))

## fit a reference model: of dd for C data only, and add to the data:--------
ref_model <- lm(data_cyano %>% filter(treat=="C"), 
                formula = pcgr ~ poly(population.mean,1)*strain + strain)
data_cyano$pcgr_ref <- predict.lm(ref_model, newdata = data_cyano)
data_cyano <- data_cyano %>%
  mutate(delta = pcgr - pcgr_ref)
#plot the data and overlay the reference model of dd
ggplot(data_cyano) +
  theme_bw() + 
  geom_point(aes(x=population.mean, y=pcgr)) + 
  geom_smooth(method=lm, aes(x=population.mean, y=pcgr), 
              formula=y ~ poly(x,1), se=F) + #
  geom_line(aes(x=population.mean, y=pcgr_ref)) +
  facet_wrap(vars(treat, strain), scales="free") #, rows = vars(ID_spec)strain ~ 

#now check if delta with ref model depends on the trait
ggplot(data_cyano) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette[c(2,3,1,4)]) + 
  geom_point(aes(x=mean, y=delta, col=as_factor(treat))) + 
  #geom_smooth(method=lm, aes(x=mean, y=delta, col=as_factor(treat)),
  #            formula=y ~ poly(x,1), se=F) + #
  facet_wrap(vars(strain), scales="free") #, rows = vars(ID_spec)strain ~ 
#It does seem to be the case. 
