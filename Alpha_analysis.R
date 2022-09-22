library(tidyverse)
library(lubridate)
library(mgcv)
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# CILIATES ------
## Import data and compute pcgr----
data_cilia <- read.csv("ciliates/DIVERCE_TdB_Ciliates_Traits.csv") %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(pcgr = log(lead(Count, 1)/Count)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Count) 

## Plot growth curves ----
ggplot(data_cilia) +
  theme_classic() + 
  aes(x=Days_fromstart, y=log10(Count), 
      col=as_factor(Temp)) + 
  geom_point() + 
  facet_grid(cols = vars(Atrazine), rows = vars(ID_spec), scales="free") #, rows = vars(ID_spec)

## Kick out: ----
#1 extremely negative pcgr: (no apparent lag phase and no real crashes)
data_cilia_clean <- data_cilia %>% 
  filter(pcgr > -5)

## Check linear growth----
ggplot(data_cilia_clean) +
  theme_classic() + 
  aes(x=Count, y=pcgr, 
      col=as_factor(Temp)) + 
  geom_point() + 
  facet_grid(cols = vars(ID_spec), rows = vars(Atrazine), scales="free") +#, rows = vars(ID_spec)
  geom_smooth(method="lm")

## Compute the alphas and analyse ----
data_cilia_clean_analysis <- data_cilia_clean %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(alpha = cov(Count, pcgr) / var(Count)) %>% #alphas
  mutate(igr = mean(pcgr) - alpha*mean(Count)) %>% #intrinsic growth rate
  mutate(cv_ar = sd_ar/mean_ar) %>% #compute cvs
  mutate(cv_area = sd_area/mean_area) %>%
  mutate(cv_speed = sd_speed/mean_speed) %>%
  mutate(cv_linearity = sd_linearity/mean_linearity) %>%
  select(-contains(c("sd_"))) %>%
  pivot_longer(contains(c("cv_", "mean_")), names_to = "trait", values_to = "value") %>%
  group_by(Atrazine, Temp, ID_spec, trait) %>%
  mutate(value_delta = max(value) - min(value)) %>% # change of mean and sd of traits
  summarise_all(mean, na.rm=T) 

ggplot(data_cilia_clean_analysis) + 
  theme_classic() + 
  aes(x=log10(value_delta), y=log10(-alpha)) + #col=as_factor(treat)
  geom_point() + 
  facet_wrap(.~trait, scales="free") + 
  geom_smooth(method="lm", se=T) #, rows = vars(ID_spec)

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
  mutate(strain_treat = paste(strain, treat, sep="_"))
#we need to make a dummy variable strain_treat because s() can only take one grouping variable

## Plot growth curves ----
ggplot(data_cyano) +
  theme_classic() + 
  aes(x=day, y=log10(population.mean), 
      col=as_factor(treat), pch=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(strain), scales="free") +
  geom_smooth(method = "gam")
#observe lag phaases and population crashes. Solution for crashes:
#fit a gam, make predictions, compute the pcgr from these predictions 
#and find the timepoint with max abundance

## Apply GAM fitting ---- 
model <- gam(log10(population.mean) ~ as_factor(strain_treat) + 
               s(day, by = as_factor(strain_treat), k=5), 
             data = data_cyano) 
data_cyano$predictions <- predict.gam(model) 
## Plot GAM predictions ----
ggplot(data_cyano) +
  theme_classic() + 
  aes(x=day, y=predictions, 
      col=as_factor(treat), pch=as_factor(treat)) + 
  geom_line() + 
  geom_point(aes(x=day, y=log10(population.mean))) + 
  facet_grid(cols = vars(strain), scales="free") #, rows = vars(ID_spec)
## Now clean based on GAM predictions: ----
# 1/remove all data for timepoints after peak abundance. 
# Note that pcgr on last day will remain there.
# 2/remove lag phase data 
data_cyano_clean <- data_cyano %>%
  group_by(strain, treat) %>%
  mutate(pcgr_pred = (lead(predictions,1) - predictions)/(lead(day,1) - day), .after = pcgr) %>%
  filter(day < max(day*(predictions==max(predictions)))) %>%
  mutate(lag_test = min(day*(ifelse(abs(pcgr_pred)>1e-1, TRUE, NA)), na.rm=T), .after = pcgr_pred) %>%
  filter(day >= lag_test)
## Plot the clean data ----
ggplot(data_cyano_clean) +
  theme_classic() + 
  aes(x=day, y=log10(population.mean), 
      col=as_factor(treat), pch=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(strain), scales="free") #, rows = vars(ID_spec)

## Now compute the alpha's and analyse----
data_cyano_analysed <- data_cyano_clean %>%
  mutate(FSC.HLin.cv = FSC.HLin.sd/FSC.HLin.mean) %>%
  mutate(SSC.HLin.cv = SSC.HLin.sd/FSC.HLin.mean) %>%
  mutate(GRN.B.HLin.cv = GRN.B.HLin.sd/GRN.B.HLin.mean) %>%
  mutate(YEL.B.HLin.cv = YEL.B.HLin.sd/YEL.B.HLin.mean) %>%
  mutate(RED.B.HLin.cv = RED.B.HLin.sd/RED.B.HLin.mean) %>%
  mutate(RED.R.HLin.cv = RED.R.HLin.sd/RED.R.HLin.mean) %>%
  select(-contains(c("NIR."))) %>% #kick out NIR data
  group_by(strain, treat) %>%
  mutate(alpha = cov(population.mean, pcgr) / var(population.mean)) %>% #alphas
  mutate(igr = mean(pcgr) - alpha*mean(population.mean)) %>% #intrinsic growth rate
  select(-contains(c(".sd", "population"))) %>%
  pivot_longer(contains(c(".cv", ".mean")), names_to = "trait", values_to = "value") %>%
  group_by(strain, treat, trait) %>%
  mutate(value_delta = max(value) - min(value)) %>% # change of mean and sd of traits
  summarise_all(mean, na.rm=T) 

ggplot(data_cyano_analysed) + 
  theme_classic() + 
  aes(x=log10(value_delta), y=log10(-alpha)) + #col=as_factor(treat)
  geom_point() + 
  facet_wrap(.~trait, scales="free") + 
  geom_smooth(method="lm", se=T) #, rows = vars(ID_spec)



