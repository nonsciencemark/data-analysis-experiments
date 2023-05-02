#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT TOOLS AND DATA -----------------------
source("Tools and data.R")
# PICK DATA SOURCE ---------------
model_system <- "cilia" #cilia or cyano
data         <- get(paste("data_",model_system, sep=""))
# DO ANALYSES --------------------
### do the stats for r and traits (t), only control --------
#data_ref <- data %>% 
#  filter(atrazine%in%c("0","no"), temperature %in%c("normal", "22")) 
stats_result <- modelling(data=data, 
                          var_to_nest_by = c("strain", "trait", "treat"),
                          formulas=c("dT ~ density + mean", #+ density*mean
                                     "dT ~ mean",
                                     "pcgr ~ mean + density",# + density*mean
                                     "pcgr ~ density")) %>%
  select(-data)

#Now join these models to all the data (just the coefficients suffice)
#and make the predictions
data_preds <- data %>%
  ungroup() %>%
  nest_by(species, strain, trait, treat) %>%
  left_join(stats_result, by=c("strain", "trait", "treat"), multiple = "all") %>%
  mutate(predictions = list(cbind(data, pred=predict.lm(model, newdata=data)))) %>%
  rowwise() %>%
  mutate(AIC=AIC(model)) %>%  
  ungroup() %>%
  mutate(response = case_when(form=="dT ~ density + mean" ~ "trait change",
                              form=="dT ~ mean" ~ "trait change",
                              form=="pcgr ~ mean + density" ~ "growth",
                              form=="pcgr ~ density" ~ "growth")) %>%
  ungroup()
  
# Plot model fit and compare (AIC)
data_preds_synth <- data_preds %>%
  mutate(predictors = case_when(form=="dT ~ density + mean" ~ "both",
                                form=="dT ~ mean" ~ "focal",
                                form=="pcgr ~ mean + density" ~ "both",
                                form=="pcgr ~ density" ~ "focal")) %>%
  select(c("species", "strain", "trait", "treat", "response", "AIC", "predictors")) %>%
  pivot_wider(names_from = "predictors", values_from = "AIC") %>%
  mutate(delta_AIC = both - focal) #if <0 then both predict better

ggplot(data_preds_synth) + 
  scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=response, y=delta_AIC, col=as.factor(treat), pch=trait) + 
  geom_jitter(width = 0.25) +
  geom_hline(yintercept = 0, lty="dotted") +
  facet_wrap(vars(species, strain), ncol=2, scales="free")
  
# Plot model fits
data_preds <- data_preds %>%
  select(-c("data", "model")) %>%
  unnest(predictions)

ggplot(data_preds %>% filter(response=="growth")) + 
  theme_bw() + 
  scale_shape_manual(values=0:10) + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=pcgr, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  facet_grid(form~trait) + #, scales="free"
  geom_abline(slope=1, intercept=0) 

ggplot(data_preds %>% filter(response=="trait change")) + 
  theme_bw() + 
  scale_shape_manual(values=0:10) + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=dT, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  #geom_smooth(method="lm", se=F)+
  facet_wrap(trait~form, scales="free") + #
  geom_abline(slope=1, intercept=0) 


  

