library('tidyverse')

source("data/cyanobacteria/generating_clean_data/adjust_traits.r")

## strains
s_pal <- setNames(c('darkblue', 'dodgerblue', 'firebrick', 'darkorange'), 
                  c( '2383', '2434', '2375', '2524'))

dat <- readRDS('data/cyanobacteria/generating_clean_data/final-pops-mon.RData') %>%
  dplyr::select(-all_of(c('date.time', 'dil', 'clusts.sp', 'clusts.st', 'spec.pops'))) %>%
  dplyr::select(-contains(c('SSC', 'var.'))) %>%
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
  rename(
    "Green-B" = "av.GRN.B",
    "NIR-B" = "av.NIR.B",
    "NIR-R" = "av.NIR.R",
    "Red-B" = "av.RED.B",
    "Red-R" = "av.RED.R",
    "Yellow-B" = "av.YEL.B",
    ) %>%
  rename(`Population density` = pop,
    Size = "av.FSC") %>%
  rowwise %>%
  mutate(
    Phycoerythrin = weighted.mean(
      c(`Green-B`, `NIR-B`, `NIR-R`, `Red-B`, `Red-R`, `Yellow-B`), weights$Phycoerythrin),
    Chlorophyll = weighted.mean(
      c(`Green-B`, `NIR-B`, `NIR-R`, `Red-B`, `Red-R`, `Yellow-B`), weights$Chlorophyll),
    Phycocyanin = weighted.mean(
      c(`Green-B`, `NIR-B`, `NIR-R`, `Red-B`, `Red-R`, `Yellow-B`), weights$Phycocyanin),
    # Phycoerythrin = av.YEL.B,
    # Chlorophyll = av.RED.B, 
    # Phycocyanin = av.RED.R
    )  %>%
  dplyr::select(-contains("-")) %>%
  mutate(strain = factor(strain, levels = names(s_pal), labels = names(s_pal))) %>% # nicer order
  arrange(strain, treat, repl, date) %>%
  as.data.frame %>%
  relocate(pcgr, .before = Size)

write.csv(dat, 'data/cyanobacteria/mono_data.csv', row.names = FALSE)
