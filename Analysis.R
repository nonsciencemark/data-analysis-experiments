#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT TOOLS AND DATA -----------------------
source("Tools and data.R")
# PICK DATA SOURCE ---------------
model_system <- "cyano" #cilia or cyano
data         <- get(paste("data_",model_system, sep=""))
# DO ANALYSES --------------------
### do the stats for r and traits (t), only control --------
data_ref <- data %>% 
  filter(atrazine%in%c("0","no"), temperature %in%c("normal", "22")) 
stats_result_r <- modelling(data=data_ref, 
                          var_to_nest_by = c("strain", "trait"),
                          formula="pcgr ~ density + mean + density*mean") %>%
  select(any_of(c("strain", "trait", "predictor", "Estimate"))) %>%
  pivot_wider(names_from = predictor, values_from = Estimate)
stats_result_t <- modelling(data=data_ref, 
                            var_to_nest_by = c("strain", "trait"),
                            formula="dT ~ density + mean + density*mean") %>%
  select(any_of(c("strain", "trait", "predictor", "Estimate"))) %>%
  pivot_wider(names_from = predictor, values_from = Estimate)
#Now join these models to all the data (just the coefficients suffice)
data_r <- data %>%
  left_join(stats_result_r, by=c("strain", "trait")) %>%
  mutate(pcgr_pred = `(Intercept)` + density.y*density.x + mean.y*mean.x + 
           `density:mean`*density.x*mean.x)

data_t <- data %>%
  left_join(stats_result_t, by=c("strain", "trait")) %>%
  mutate(dT_pred = `(Intercept)` + density.y*density.x + mean.y*mean.x + 
           `density:mean`*density.x*mean.x)

ggplot(data_r) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=pcgr, y= pcgr_pred, col=atrazine) + 
  geom_point() + 
  facet_grid(temperature~strain, scales="free") + 
  geom_abline(slope=1, intercept=0) 

ggplot(data_t) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=dT, y= dT_pred, col=atrazine, pch=trait) + 
  geom_point() + 
  facet_grid(temperature~strain, scales="free") +
  geom_abline(slope=1, intercept=0)

