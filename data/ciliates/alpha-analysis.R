# packages
library('tidyverse')  

# load data
data <- read.csv('ciliates/DIVERCE_TdB_Ciliates_Traits.csv') %>%
  group_by(ID_spec, Temp, Atrazine) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  select(contains(c('sd_', 'ID', 'Temp', 'Atrazine', 'r1', 'a11'))) %>%
  pivot_longer(contains(c('sd_')), names_to = 'trait', values_to = 'sd') %>%
  dplyr::filter(!grepl('Spiro', ID_spec))

# plot
ggplot(data %>% 
         mutate(Temp = as.factor(Temp), Atrazine = as.factor(Atrazine))) + 
  theme_bw() + 
  aes(x = sd, y = -a11, col = Temp, pch = Atrazine) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point() + 
  labs(x = 'Trait standard deviation', 
       y = expression(paste('Self-limitation (', alpha['11'], ')'))) +
  facet_wrap(.~trait, scales='free') #, rows = vars(ID_spec)

