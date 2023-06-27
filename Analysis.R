# IMPORT TOOLS AND DATA ----
source("Tools and data.R")
source("Plotting functions.R")
library("lmtest") # for granger causality

# PLOT BASIC OUTPUTS ----
basic_plots("cyano")
basic_plots("cilia")

## Results ----

# Ciliates: 
# negative dd, but different for different treatments; 
# negative td, not different for different treatments.

# Cyanos:
# dd (mostly negative); different for different treatments. 
# td (mostly negative); different for different treatments.

# Task 1: Predict next time point within condition: single vs. two variables ----

LHSs <- c('dT', 'pcgr')
RHSs <- c('density', 
          'trait', 
          # 'poly(trait, 2)', 
          'density + trait' # , 
          # 'density + poly(trait, 2)'
          )

eqs <- expand.grid('LHS' = LHSs, 'RHS' = RHSs) %>%
  # ensure that we're not having e.g. pcgr ~ trait
  dplyr::filter(LHS == 'dT' & grepl('trait', RHS) | 
                  LHS == 'pcgr' & grepl('density', RHS)) %>%
  apply(1, function(i) {paste0(i[1], ' ~ ', i[2])}) %>%
  sort()

AIC_plots("cyano")
AIC_plots("cilia")

granger_plots("cyano")
granger_plots("cilia")

## Results ----

# Ciliates: 
# Combining two predictors rarely helped to predict growth or trait change

# Cyanos:
# Combining two predictors rarely helped to predict growth, but did improve trait change prediction

## Task 2: Predict next time point across conditions: single vs. two predictors-----

third_plots("cyano")
third_plots("cilia")

## Results ----

# Ciliates: 
# Combining two predictors rarely helped to predict growth or trait change to other conditions

# Cyanos:
# Combining two predictors rarely helped to predict growth but could improve or worsen trait change prediction.

# Conclusion: not a very strong link between traits and abundance, with some exceptions for the cyanos

## Basic stats: treatment effects on dd of pcgr and td of dT -----

td_general_plots("cyano")
td_general_plots("cilia")
