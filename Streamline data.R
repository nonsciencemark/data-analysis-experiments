
#Import and make uniform the ciliate data
data_cilia <- read.csv("ciliates/OverviewTraitsFull.csv") %>%
  rename(density = Parts_permL) %>%
  rename(strain = ID_spec) %>%
  rename(temperature = Temp) %>%
  rename(atrazine = Atrazine) %>%
  group_by(atrazine, temperature, strain) %>%
  mutate(pcgr = log(lead(density, 1)/density)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = density) %>% 
  filter(pcgr > -5, strain!="Spiro_5") %>% #Kick out extremely negative pcgr, as well as Spiro_5 (not enough data to make a ref model)
  ungroup() %>%
  separate(strain, into = c("species", "sp.strain"), remove=F) %>%#link to species
  group_by(strain) %>%
  mutate(phase = cut_interval(log10(density), n=5, labels = FALSE)) %>%
  pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  separate(trait, into=c("stat", "trait"), sep="_") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean)

#Import and make uniform the cyano data
data_cyano <- read.csv("cyanobacteria/mono_data.csv") %>%
  rename(density = population.mean) %>%
  separate(date.time, sep=" ", into = c("date", "time")) %>%
  mutate(day = yday(date), .after = date) %>%
  select(-c(date, time)) %>%
  mutate(species = case_when(strain%in%c(2375,2524)~"one",
                             TRUE ~ "two")) %>%#link to species
  mutate(atrazine = case_when(treat%in%c("A", "AT")~"yes",
                              TRUE ~ "no")) %>%
  mutate(temperature = case_when(treat%in%c("A", "C")~"normal",
                          TRUE ~ "hot")) %>%
  group_by(atrazine, temperature, strain) %>%
  mutate(day = day - min(day)) %>%
  mutate(pcgr = log(lead(density,1)/density)/(lead(day,1)-day),
         .after = density) %>%
  filter(!is.na(pcgr)) %>%
  ungroup() %>%
  group_by(strain) %>%
  mutate(phase = cut_interval(log10(density), n=5, labels = FALSE)) %>%
  select(-population.sd) %>%
  pivot_longer(FSC.HLin.mean:NIR.R.HLin.sd, names_to="trait") %>%
  separate(trait, into=c("trait", "stat"), sep="in.") %>%
  pivot_wider(names_from=stat, values_from=value) %>%
  mutate(cv=sd/mean) 

