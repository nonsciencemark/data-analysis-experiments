library('tidyverse')

## strains
s_pal <- setNames(c('darkblue', 'dodgerblue', 'firebrick', 'darkorange'), 
                  c( '2383', '2434', '2375', '2524'))

dat <- readRDS('data/cyanobacteria/generating_clean_data/final-pops-mon.RData') %>%
  dplyr::select(-all_of(c('date.time', 'dil', 'clusts.sp', 'clusts.st', 'spec.pops'))) %>%
  dplyr::select(-contains('var.')) %>%
  dplyr::select(-contains(c('SSC', 'GRN', 'NIR'))) %>%
  dplyr::filter(date != 0, date < 11) %>%
  mutate(date = as.numeric(date + 0)) %>%
  group_by(treat, strain, repl) %>%
  mutate(pcgr = log(lead(pop, 1) / pop) / as.numeric(lead(date, 1) - date), .after = pcgr) %>%
  ungroup() %>%
  group_by(strain, treat) %>%
  mutate(pcgr_clean = ifelse(
    (pcgr < (quantile(pcgr, .25, na.rm = TRUE) - 1.5 * IQR(pcgr, na.rm = TRUE)) | 
       pcgr > (quantile(pcgr, .75, na.rm = TRUE) + 1.5 * IQR(pcgr, na.rm = TRUE))), 
    NA, pcgr), .after = pcgr) %>%
  dplyr::select(-pcgr) %>%
  rename(pcgr = pcgr_clean) %>%
  ungroup() %>%
  rename(`Population density` = pop, 
         Size = av.FSC, 
         Phycoerythrin = av.YEL.B, 
         Chlorophyll = av.RED.B, 
         Phycocyanin = av.RED.R)  %>%
  mutate(strain = factor(strain, levels = names(s_pal), labels = names(s_pal))) %>% # nicer order
  arrange(strain, treat, repl, date) %>%
  as.data.frame

write.csv(dat, 'data/cyanobacteria/mono_data.csv', row.names = FALSE)
