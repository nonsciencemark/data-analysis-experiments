# packages
library('tidyverse')  

# load data
data <- read.csv('ciliates/DIVERCE_TdB_Ciliates_Traits.csv') %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  select(contains(c('sd_', 'ID', 'Temp', 'Atrazine', 'r1', 'a11'))) %>%
  pivot_longer(contains(c('sd_')), names_to = 'trait', values_to = 'sd') %>%
  dplyr::filter(!grepl('Spiro', ID_spec)) %>%
  mutate(trait = replace(trait, trait == 'sd_ar', 'Aspect ratio'),
         trait = replace(trait, trait == 'sd_area', 'Area'),
         trait = replace(trait, trait == 'sd_speed', 'Speed'),
         trait = replace(trait, trait == 'sd_linearity', 'Linearity'))

# plot with just the traits as facets
ggplot(data %>% 
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_smooth(method = 'lm') +
  geom_point(aes(col = Temp, pch = Atrazine)) + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Self-limitation (', alpha['11'], ')'))) +
  facet_wrap(.~trait, scales = 'free')

# plot including species as facet
ggplot(data %>% 
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11, col = Temp, pch = Atrazine) + 
  scale_x_log10() + 
  scale_y_log10(breaks = trans_breaks('log10', function(x){10^x}),
                labels = trans_format('log10', scales::math_format(10^.x))) + 
  geom_point() + 
  labs(x = 'Trait standard deviation (SD)', 
       y = expression(paste('Self-limitation (', alpha['11'], ')'))) +
  facet_grid(ID_spec~trait)
