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
data_cilia <- read.csv("data/ciliates/DIVERCE_TdB_Ciliates_Traits_FULL.csv") %>%
  rename(density = Parts_permL) %>%
  rename(treat = Treatment) %>%
  mutate(strain = as.factor(ID_spec)) %>%
  mutate(temperature = as.factor(Temp)) %>%
  mutate(atrazine = as.factor(Atrazine)) %>%
  group_by(atrazine, temperature, strain) %>%
  mutate(pcgr = log(lead(density, 1)/density)/(lead(Time_Days, 1)-Time_Days),
         .after = density) %>% 
  # filter(!is.na(pcgr)) %>%
  filter(pcgr > -5, strain!="Spiro_5") %>% #Kick out extremely negative pcgr, as well as Spiro_5 (not enough data to make a ref model)
  ungroup() %>%
  separate(strain, into = c("species", "sp.strain"), remove=F) #link to species

data_cilia_pca <- data_cilia %>%
  group_by(strain) %>%
  nest() %>%
  mutate(pca_data = map(data, ~.x %>% dplyr::select(contains("mean")))) %>% 
  mutate(pca_scores = map2(data, pca_data, ~cbind(.x, pca1=prcomp(.y, center=T, scale=T)$x[,1]))) %>%
  dplyr::select(c("strain", "pca_scores")) %>%
  unnest(cols = c(pca_scores))

data_cilia <- data_cilia_pca %>%
  rename(trait=pca1) %>%
  group_by(strain, atrazine, temperature) %>%
  mutate(dT = lead(trait, 1)-trait/(lead(Time_Days, 1)-Time_Days)) %>%
  filter(!is.na(dT)) %>%
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
data_cilia$strain <- factor(data_cilia$strain, 
                            levels=c("Spiro_C", "Spiro_D", "Tetra_1", "Tetra_2", 
                                     "Loxo_1", "Loxo_2", "Para_4"))

#Import and make uniform the cyano data-------------------
data_cyano <- read.csv("data/cyanobacteria/mono_data.csv") %>%
  rename(density = Population.density, day = date) %>%
  # separate(date.time, sep=" ", into = c("date", "time")) %>%
  # mutate(day = yday(date), .after = date) %>%
  # dplyr::select(-c(date, time)) %>%
  # fixed species names
  mutate(species = case_when(strain %in% c(2375, 2524) ~ "V", TRUE ~ "VIII")) %>% #link to species
  mutate(strain = as_factor(strain)) %>%
  mutate(atrazine = as_factor(case_when(treat %in% c("C", "T") ~ "no", TRUE ~ "yes"))) %>%
  mutate(temperature = as_factor(case_when(treat %in% c("A", "C") ~ "normal", TRUE ~ "hot"))) %>%
  group_by(atrazine, temperature, strain) %>%
  # mutate(day = day - min(day)) %>%
  # mutate(pcgr = log(lead(density,1)/density)/(lead(day,1)-day),
  #        .after = density) %>%
  # filter(!is.na(pcgr)) %>% # removed this b/c it gets done later and it was causing data loss
  ungroup() #%>%
  # dplyr::select(-population.sd) 

data_cyano_pca <- data_cyano %>%
  group_by(strain) %>%
  # not really a good idea to include non-functional traits
  # dplyr::select(-c('SSC.HLin.mean', 'NIR.B.HLin.mean', 'NIR.R.HLin.mean', 'GRN.B.HLin.mean')) %>%
  nest() %>%
  mutate(pca_data = map(data, ~.x %>% dplyr::select(Size:Phycocyanin))) %>% 
  mutate(pca_scores = map2(data, pca_data, ~cbind(.x, pca1=prcomp(.y, center=T, scale=T)$x[,1]))) %>%
  dplyr::select(c("strain", "pca_scores")) %>%
  unnest(cols = c(pca_scores))
  
data_cyano <- data_cyano_pca %>%
  # fixed here too
  # dplyr::select(-c("FSC.HLin.mean":"YEL.B.HLin.sd")) %>%
  rename(trait=pca1) %>%
  #group_by(strain, trait, atrazine, temperature) %>%
  group_by(strain, atrazine, temperature, repl) %>%
  #mutate(dT = lead(mean, 1)-mean/(lead(day, 1)-day)) %>%
  mutate(dT = lead(trait, 1)-trait/(lead(day, 1)-day)) %>%
  filter(!is.na(dT)) %>% 
  mutate(treat = as.factor(treat)) 
  
data_cyano$treat <- factor(data_cyano$treat, levels=c("C", "T", "A", "AT"))
# this is more consistent with ciliate method
data_cyano$strain <- paste(data_cyano$species, "_", data_cyano$strain, sep="")

