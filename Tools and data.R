library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
#Functions -----------------
modelling <- function(data=data, var_to_nest_by="strain", formulas) {
  var_to_nest_by <- c(var_to_nest_by, "form")
  test <- expand_grid(data%>%ungroup(), 
                      form = formulas) %>%
    nest_by(across({{var_to_nest_by}})) %>%
    mutate(model = list(lm(as.formula(form), 
                                   data = data)))
  return(test)
}

#Import and make uniform the ciliate data ----------
data_cilia <- read.csv("data/ciliates/OverviewTraitsFull.csv") %>%
  rename(density = Parts_permL) %>%
  rename(treat = Treatment) %>%
  mutate(strain = as.factor(ID_spec)) %>%
  mutate(temperature = as.factor(Temp)) %>%
  mutate(atrazine = as.factor(Atrazine)) %>%
  group_by(atrazine, temperature, strain) %>%
  mutate(pcgr = log(lead(density, 1)/density)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = density) %>% 
  filter(pcgr > -5, strain!="Spiro_5") %>% #Kick out extremely negative pcgr, as well as Spiro_5 (not enough data to make a ref model)
  ungroup() %>%
  separate(strain, into = c("species", "sp.strain"), remove=F) %>%#link to species
  pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  separate(trait, into=c("stat", "trait"), sep="_") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean) %>%
  group_by(strain, trait, atrazine, temperature) %>%
  mutate(dT = lead(mean, 1)-mean/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = cv) %>%
  filter(Temp>20) %>%
  mutate(treat = case_when((Atrazine==0)&(Temp==22) ~ "C",
                           (Atrazine==10)&(Temp==22) ~ "a",
                           (Atrazine==20)&(Temp==22) ~ "A",
                           (Atrazine==0)&(Temp==24) ~ "T",
                           (Atrazine==10)&(Temp==24) ~ "aT",
                           (Atrazine==20)&(Temp==24) ~ "AT"))

data_cilia$species <- factor(data_cilia$species, levels=c("Loxo", "Spiro", "Tetra", "Para"))
data_cilia$treat <- factor(data_cilia$treat, levels=c("C", "T", 
                                                      "a", "A",
                                                      "aT", "AT"))

#Import and make uniform the cyano data-------------------
data_cyano <- read.csv("data/cyanobacteria/mono_data.csv") %>%
  rename(density = population.mean) %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  mutate(species = case_when(strain%in%c(2375,2524)~"one",
                             TRUE ~ "two")) %>%#link to species
  mutate(strain = as_factor(strain)) %>%
  mutate(atrazine = as_factor(case_when(treat%in%c("C", "T")~"no",
                              TRUE ~ "yes"))) %>%
  mutate(temperature = as_factor(case_when(treat%in%c("A", "C")~"normal",
                          TRUE ~ "hot"))) %>%
  group_by(atrazine, temperature, strain) %>%
  mutate(day = day - min(day)) %>%
  mutate(pcgr = log(lead(density,1)/density)/(lead(day,1)-day),
         .after = density) %>%
  filter(!is.na(pcgr)) %>%
  ungroup() %>%
  select(-population.sd) %>%
  pivot_longer(FSC.HLin.mean:NIR.R.HLin.sd, names_to="trait") %>%
  separate(trait, into=c("trait", "stat"), sep="in.") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean) %>% 
  group_by(strain, trait, atrazine, temperature) %>%
  mutate(dT = lead(mean, 1)-mean/(lead(day, 1)-day),
         .after = cv) %>%
  mutate(treat = as.factor(treat)) 
  
data_cyano$treat <- factor(data_cyano$treat, levels=c("C", "T", "A", "AT"))


