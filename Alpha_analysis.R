library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)
# changed by TdB for trait index, do not quite get what I did yet...

# CILIATES ------
## Import data and compute pcgr, and only keep specific species-trait combinations----
data_cilia <- read.csv("C:/Users/debruin/OneDrive - UCL/Bio/UCL/GitHub/TeamWork/data-analysis-experiments/ciliates/OverviewTraitsFull.csv") %>%
  mutate(strain = as_factor(ID_spec)) %>%
  group_by(Atrazine, Temp, strain) %>%
  mutate(pcgr = log(lead(Parts_permL, 1)/Parts_permL)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Parts_permL) %>% 
  filter(pcgr > -5, strain!="Spiro_5") %>% #Kick out extremely negative pcgr, as well as Spiro_5 (not enough data to make a ref model)
  ungroup() %>%
  separate(strain, into = c("species", "sp.strain"), remove=F) %>%#link to species
  group_by(strain) %>%
  mutate_at(c("mean_speed", "mean_area", "mean_ar", "mean_linearity", "sd_speed", "sd_area", "sd_ar", "sd_linearity", "pcgr", "Parts_permL", "Days_fromstart"), ~ (scale(.) %>% as.vector)) %>%
  mutate(phase = cut_interval(Parts_permL, n=5, labels = FALSE)) %>%
  ungroup() %>%
  mutate(Strategy_Go = (mean_speed + mean_linearity), 
         Strategy_Stay = (mean_area + mean_ar),
         Strategy_index = (-mean_speed + -mean_linearity + mean_area + mean_ar)) %>% # index by summing speed, linearity, size and shape
  pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  separate(trait, into=c("stat", "trait"), sep="_") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean)
#%>%
#  filter(((species=="Spiro")&(trait=="speed"))|((species=="Tetra")&(trait=="speed"))
#         |((species=="Para")&(trait=="speed"))|((species=="Loxo")&(trait=="speed")))

##check if dd changes with treatment -------
ggplot(data_cilia) +
  scale_colour_manual(values=c("black", cbPalette)) + 
  aes(x=Parts_permL,y=pcgr, col=as.factor(Treatment)) +
  geom_point() +
  facet_grid(cols=vars(strain), scales="free") +
  geom_smooth(method=lm, aes(x=Parts_permL, y=pcgr, col=as.factor(Treatment)),
              formula=y ~ poly(x,1), se=F) 

##check if dependence of trait on density changes with treatment -------
ggplot(data_cilia) +
  scale_colour_manual(values=c("black", cbPalette)) + 
  aes(x=(Parts_permL),y=(Strategy_index), col=as.factor(Atrazine)) +
  geom_point() +
  facet_grid(cols=vars(strain), rows=vars(Temp), scales="free") +
  geom_smooth(method=lm, aes(x=(Parts_permL), y=(Strategy_index), col=as.factor(Atrazine)),
              formula=y ~ poly(x,1), se=F)

## fit a reference model of dd and trait dependence on density--------
# This model only uses the control data
# w/o crashes. Crash detection happens with a GAM
# No apparent lag phases here.
data_cilia_ref <- data_cilia %>% 
  filter(Treatment==1) 
gam_model <- gam(Parts_permL ~ strain + 
                   s(Days_fromstart, by = strain, k=5),
                 data = data_cilia_ref)
data_cilia_ref$predictions <- predict.gam(gam_model) 
data_cilia_ref <- data_cilia_ref %>%
  group_by(strain) %>%
  filter(Days_fromstart < max(Days_fromstart*(predictions==max(predictions))))
###first dd----------------
ref_model_dd <- lm(data_cilia_ref, formula = pcgr ~ poly(Parts_permL,1)*strain + strain)

###then traits----------------
ref_model_trait <- lm(data_cilia_ref,  
                      formula = mean ~ poly(Parts_permL,1)*strain*Strategy_index + strain*Strategy_index + strain + Strategy_index)
###add predictions to the original data frame------------
data_cilia$pcgr_ref <- predict.lm(ref_model_dd, newdata = data_cilia)
data_cilia$mean_ref <- predict.lm(ref_model_trait, newdata = data_cilia)
data_cilia <- data_cilia %>%
  mutate(delta_pcgr = pcgr-pcgr_ref) %>%
  mutate(delta_index = mean-mean_ref)

##now check if delta with ref model depends on the trait----------
ggplot(data_cilia) +
  theme_bw() + 
  scale_size_manual(values=c(rep(1,2), 3, rep(1,10))) +
  scale_colour_manual(values=c("black", cbPalette)) + 
  geom_point(aes(x=delta_index, y=delta_pcgr, col=as.factor(Atrazine)), 
             alpha=0.8) + 
  geom_smooth(method=lm, aes(x=delta_index, y=delta_pcgr, col=as.factor(Atrazine)),#
              formula=y ~ poly(x,1)) + #,  col="black"
  facet_grid(rows=vars(species), cols=vars(Temp), scales="free")

# CYANO ------
## Import data and compute pcgr----
data_cyano <- read.csv("C:/Users/debruin/OneDrive - UCL/Bio/UCL/GitHub/TeamWork/data-analysis-experiments/cyanobacteria/mono_data.csv") %>%
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
  mutate(cv=sd/mean) 
#%>%
#filter(((species=="one")&(trait=="RED.B.HL"))|((species=="two")&(trait=="YEL.B.HL")))

##check if dd changes with treatment -------
ggplot(data_cyano) +
  aes(x=population.mean,y=pcgr, col=treat) +
  geom_point() +
  facet_grid(cols=vars(strain), scales="free") +
  geom_smooth(method=lm, aes(x=population.mean, y=pcgr, col=treat),
              formula=y ~ poly(x,1), se=F) 

##check if dependence of trait on density changes with treatment -------
ggplot(data_cyano) +
  aes(x=log10(population.mean),y=mean, col=treat) +
  geom_point() +
  facet_grid(cols=vars(strain), rows=vars(trait), scales="free") +
  geom_smooth(method=lm, aes(x=log10(population.mean), y=mean, col=treat),
              formula=y ~ poly(x,2), se=F)

#yes, 2nd order

## fit a reference model of dd and trait dependence on density--------
# This model only uses the control data
# no apparent lags; no apparent crashes
###first dd----------------
ref_model_dd <- lm(data_cyano%>%filter(treat=="C"), formula = pcgr ~ poly(population.mean,1)*strain + strain)

###then traits----------------
ref_model_trait <- lm(data_cyano%>%filter(treat=="C"),  
                      formula = mean ~ poly(log10(population.mean),2)*strain*trait + strain*trait + strain + trait)

###add predictions to the original data frame -------------
data_cyano$pcgr_ref <- predict.lm(ref_model_dd, newdata = data_cyano)
data_cyano$mean_ref <- predict.lm(ref_model_trait, newdata = data_cyano)
data_cyano <- data_cyano %>%
  mutate(delta_pcgr = pcgr-pcgr_ref) %>%
  mutate(delta_trait = log10(mean/mean_ref))

##now check if delta with ref model depends on the trait----------
ggplot(data_cyano) +
  theme_bw() + 
  scale_size_manual(values=c(rep(1,2), 3, rep(1,10))) +
  scale_colour_manual(values=cbPalette[c(2,3,1,4:8)]) + 
  geom_point(aes(x=delta_trait, y=delta_pcgr, col=treat, pch=strain), 
             alpha=0.8) + 
  geom_smooth(method=lm, aes(x=delta_trait, y=delta_pcgr, col=treat),#
              formula=y ~ poly(x,1)) + #,  col="black"
  facet_grid(rows=vars(species), cols=vars(trait), scales="free")


