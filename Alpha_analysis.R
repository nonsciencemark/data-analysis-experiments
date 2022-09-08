library(tidyverse)
library(lubridate)

# CILIATES ------
## Import data and compute pcgr----
data_cilia <- read.csv("ciliates/DIVERCE_TdB_Ciliates_Traits.csv") %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(Days_fromstart.0 = lag(Days_fromstart, 1), .after = Days_fromstart) %>%
  mutate(Count.0 = lag(Count, 1), 
         .after = Count) %>%
  mutate(pcgr = log(Count/Count.0)/(Days_fromstart-Days_fromstart.0),
         .after = Count.0) 

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
  filter(pcgr > -10)

## Check linear growth----
ggplot(data_cilia_clean) +
  theme_classic() + 
  aes(x=Count.0, y=pcgr, 
      col=as_factor(Temp)) + 
  geom_point() + 
  facet_grid(cols = vars(ID_spec), rows = vars(Atrazine), scales="free") +#, rows = vars(ID_spec)
  geom_smooth(method="lm")

## Continue with the analysis----
data_cilia_clean_analysis <- data_cilia_clean %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(alpha = cov(Count.0, pcgr) / var(Count.0)) %>% #alphas
  mutate(cv_ar = sd_ar/mean_ar) %>% #compute cvs
  mutate(cv_area = sd_area/mean_area) %>%
  mutate(cv_speed = sd_speed/mean_speed) %>%
  mutate(cv_linearity = sd_linearity/mean_linearity) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("cv_", "ID", "Temp", "Atrazine", "alpha"))) %>%
  pivot_longer(contains(c("cv_")), names_to = "trait", values_to= "cv") 

ggplot(data_cilia_clean_analysis) + 
  theme_classic() + 
  aes(x=log10(cv), y=log(-alpha), 
      col=as_factor(Temp), pch=as_factor(Atrazine)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") + #, rows = vars(ID_spec)
  geom_smooth(method="lm", se=F)

# CYANO ------
## Import data and compute pcgr----
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
         .after = population.mean.0) %>%
  mutate(strain_treat = paste(strain, treat, sep="_"))
#we need to make a dummy variable because s() can only take one grouping variable

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
## Now clean based on GAM predictions:
# remove all data for timepoints after peak abundance
# 
data_cyano_clean <- data_cyano %>%
  group_by(strain, treat) %>%
  filter(day < max(day*(predictions==max(predictions))))
  mutate(lag.phase = ((log10(population.mean.0)<3.5)&(log10(population.mean)<3.5))) %>%
  filter(lag.phase == FALSE)

#all data points where log10(population.mean) is consistently low (<3.5) because of lags


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
  mutate(GRN.B.HLin.cv = GRN.B.HLin.sd/GRN.B.HLin.mean) %>%
  mutate(YEL.B.HLin.cv = YEL.B.HLin.sd/YEL.B.HLin.mean) %>%
  mutate(RED.B.HLin.cv = RED.B.HLin.sd/RED.B.HLin.mean) %>%
  mutate(RED.R.HLin.cv = RED.R.HLin.sd/RED.R.HLin.mean) %>%
  group_by(strain, treat) %>%
  mutate(alpha = cov(population.mean.0, pcgr) / var(population.mean.0)) %>% #alphas
  group_by(strain, treat) %>%
  summarise_all(mean, na.rm=T) %>%
  select(contains(c("strain", "treat", "alpha", ".cv"))) %>%
  pivot_longer(contains(c(".cv")), names_to = "trait", values_to= "cv")

ggplot(data_cyano %>% filter(treat %in% c("C", "T"))) + 
  theme_classic() + 
  aes(x=log10(cv), y=alpha, 
      col=as_factor(treat)) + 
  geom_point() + 
  facet_grid(cols = vars(trait), scales="free") + 
  geom_smooth(method="lm", se = F) #, rows = vars(ID_spec)



