library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# CILIATES ------
## Import data and compute pcgr, and only keep specific species-trait combinations----
data_cilia <- read.csv("ciliates/OverviewTraitsFull.csv") %>%
  rename(strain = ID_spec) %>%
  group_by(Atrazine, Temp, strain) %>%
  mutate(pcgr = log(lead(Parts_permL, 1)/Parts_permL)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Parts_permL) %>% 
  filter(pcgr > -5, strain!="Spiro_5") %>% #Kick out extremely negative pcgr, as well as Spiro_5 (not enough data to make a ref model)
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

## fit a reference model for density dependence ------
# This model only uses the control data, cleaned up, i.e. 
# w/o crashes and lag phases if any. Crash detection happens with a GAM
# No apparent lag phases here.
data_cilia_ref <- data_cilia %>% 
  filter(Treatment==1) 

gam_model <- gam(log10(Parts_permL) ~ as_factor(strain) + 
             s(Days_fromstart, by = as_factor(strain), k=5),
             data = data_cilia_ref)
data_cilia_ref$predictions <- predict.gam(gam_model) 
data_cilia_ref <- data_cilia_ref %>%
  group_by(strain) %>%
  filter(Days_fromstart < max(Days_fromstart*(predictions==max(predictions))))

ref_model <- lm(data_cilia_ref, 
                formula = pcgr ~ poly(Parts_permL,1)*strain + strain)
#add to the original data frame
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
ggplot(data_cilia) +
  theme_bw() + 
  scale_size_manual(values=c(3, rep(1,10))) +
  scale_colour_manual(values=c(cbPalette, "black")) + 
  geom_point(aes(x=mean, y=delta, col=as_factor(Treatment),
                 size=as_factor(Treatment)),
             alpha=0.5) + 
  #geom_smooth(method=lm, aes(x=mean, y=delta, col=as_factor(Treatment)),
  #            formula=y ~ poly(x,2), se=F) + #
  facet_wrap(vars(strain), scales="free") #, rows = vars(ID_spec)strain ~ 
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

## fit a reference model of dd--------
# This model only uses the control data, cleaned up, i.e. 
# w/o crashes and lag phases if any. Crash detection happens with a GAM
# no apparent lags
data_cyano_ref <- data_cyano %>% 
  filter(treat=="C") 

gam_model <- gam(log10(population.mean) ~ as_factor(strain) + 
                   s(day, by = as_factor(strain), k=5),
                 data = data_cyano_ref)
data_cyano_ref$predictions <- predict.gam(gam_model) 
data_cyano_ref <- data_cyano_ref %>%
  group_by(strain) %>%
  filter(day < max(day*(predictions==max(predictions))))

ref_model <- lm(data_cyano_ref, formula = pcgr ~ poly(population.mean,1)*strain + strain)
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
  scale_size_manual(values=c(rep(1,2), 3, rep(1,10))) +
  scale_colour_manual(values=cbPalette[c(2,3,1,4)]) + 
  geom_point(aes(x=mean, y=delta, col=as_factor(treat), 
                 size=as_factor(treat)), 
             alpha=0.8) + 
  #geom_smooth(method=lm, aes(x=mean, y=delta, col=as_factor(treat)),
  #            formula=y ~ poly(x,1), se=F) + #
  facet_wrap(vars(strain), scales="free") #, rows = vars(ID_spec)strain ~ 
#It does seem to be the case. 

# To do:
# Standardize traits
# Do clustering algorithm cf. Lisa's paper? 

