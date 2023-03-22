
#Functions -----------------
modelling <- function(data=data, var_to_nest_by="strain", formula) {
  test <- data %>%
    ungroup()%>%
    nest_by(across({{var_to_nest_by}})) %>%
    mutate(model = list(lm(as.formula(formula), 
                                   data = data))) %>%
    mutate(model_summary = 
             list(as_tibble(rownames_to_column(as.data.frame(summary(model)$coefficients),
                                               var = "predictor")))) %>%
    unnest(model_summary) %>%
    select(-c("data")) 
  return(test)
}

#Import and make uniform the ciliate data ----------
data_cilia <- read.csv("data/ciliates/OverviewTraitsFull.csv") %>%
  rename(density = Parts_permL) %>%
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
  filter(Temp>20)

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
         .after = cv)
  
data_cyano$atrazine <- relevel(data_cyano$atrazine, ref="no")

