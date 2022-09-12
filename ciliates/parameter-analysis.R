# packages
library('tidyverse')  
<<<<<<< Updated upstream
=======
library('ggplot2')
>>>>>>> Stashed changes
library('scales')

# load data ====
# data <- read.csv('ciliates/DIVERCE_TdB_Ciliates_Traits.csv') %>%
#   group_by(ID_spec, Temp, Atrazine) %>%
#   summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
#   select(contains(c('sd_', 'ID', 'Temp', 'Atrazine', 'r1', 'a11'))) %>%
#   pivot_longer(contains(c('sd_')), names_to = 'trait', values_to = 'sd') %>%
#   dplyr::filter(!grepl('Spiro', ID_spec)) %>%
#   mutate(trait = replace(trait, trait == 'sd_ar', 'Aspect ratio'),
#          trait = replace(trait, trait == 'sd_area', 'Area'),
#          trait = replace(trait, trait == 'sd_speed', 'Speed'),
#          trait = replace(trait, trait == 'sd_linearity', 'Linearity'))

data1 <- read.csv('C:/Users/debruin/OneDrive - UCL/Bio/UCL/GitHub/TeamWork/data-analysis-experiments/ciliates/DIVERCE_TdB_Ciliates_Traits.csv') %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  select(contains(c('sd_', 'mean_', 'ID', 'Temp', 'Atrazine', 'r1', 'a11'))) %>%
  pivot_longer(contains(c('sd_', 'mean_')), 
               names_to = c('metric', 'trait'), 
               names_pattern = '(.*)_(.*)') %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  #dplyr::filter(!grepl('Spiro', ID_spec)) %>%
  mutate(trait = replace(trait, trait == 'sd_ar', 'Aspect ratio'),
         trait = replace(trait, trait == 'sd_area', 'Area'),
         trait = replace(trait, trait == 'sd_speed', 'Speed'),
         trait = replace(trait, trait == 'sd_linearity', 'Linearity'))

data1 <- data1[-c(1:12),]

data2 <- read.csv('C:/Users/debruin/OneDrive - UCL/Bio/UCL/GitHub/TeamWork/data-analysis-experiments/ciliates/DIVERCE_TdB_Ciliates_Traits_AR_filt_Loxo.csv') %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  select(contains(c('sd_', 'mean_', 'ID', 'Temp', 'Atrazine', 'r1', 'a11'))) %>%
  pivot_longer(contains(c('sd_', 'mean_')), 
               names_to = c('metric', 'trait'), 
               names_pattern = '(.*)_(.*)') %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  #dplyr::filter(!grepl('Spiro', ID_spec)) %>%
  mutate(trait = replace(trait, trait == 'sd_ar', 'Aspect ratio'),
         trait = replace(trait, trait == 'sd_area', 'Area'),
         trait = replace(trait, trait == 'sd_speed', 'Speed'),
         trait = replace(trait, trait == 'sd_linearity', 'Linearity'))

# alpha analysis ====
dataOUT <- data %>%
  filter(a11 < -1e-8)
# plot aii with just the traits as facets
ggplot(data1 %>% 
         #filter(ID_spec == "Loxo_1") %>%
         mutate(Temp = as.factor(ID_spec), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_smooth(method = 'lm') +
  geom_point(aes(col = Temp, pch = Atrazine)) + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Self-limitation (', alpha['ii'], ')'))) +
  facet_wrap(.~trait, scales = 'free')

ggplot(data2 %>% 
         #filter(ID_spec == "Loxo_1") %>%
         mutate(Temp = as.factor(ID_spec), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_smooth(method = 'lm') +
  geom_point(aes(col = Temp, pch = Atrazine)) + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Self-limitation (', alpha['ii'], ')'))) +
  facet_wrap(.~trait, scales = 'free')

# plot aii including species as facet
ggplot(dataOUT %>% 
         filter(ID_spec == "Loxo_1" | ID_spec == "Loxo_2") %>%
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11, col = Temp, pch = Atrazine) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_point() + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Self-limitation (', alpha['ii'], ')'))) +
  facet_grid(ID_spec~trait)

# mu analysis ==== 

# plot ri with just the traits as facets
ggplot(dataOUT %>% 
         filter(ID_spec == "Loxo_1" | ID_spec == "Loxo_2") %>%
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = mean, y = r1) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_smooth(method = 'lm') +
  geom_point(aes(col = Temp, pch = Atrazine)) + 
  labs(x = 'Trait mean', 
       y = expression(paste('Intrinsic growth rate (', mu['i'], ')'))) +
  facet_wrap(.~trait, scales = 'free')

# plot ri including species as facet
ggplot(dataOUT %>% 
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = mean, y = r1, col = Temp, pch = Atrazine) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_point() + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Intrinsic growth rate (', mu['i'], ')'))) +
  facet_grid(ID_spec~trait)

