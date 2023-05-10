#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT TOOLS AND DATA -----------------------
source("Tools and data.R")
# PICK DATA SOURCE ---------------
model_system <- "cyano" #cilia or cyano
data         <- get(paste("data_",model_system, sep=""))
# DO ANALYSES --------------------
## Check correlations between traits and abundance -----
ggplot(data) + 
  scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=log10(density), y=trait, col=as.factor(treat)) + 
  #aes(x=response, y=delta_AIC, col=as.factor(treat), pch=trait) + 
  #geom_jitter(width = 0.25) +
  geom_point(pch=1) +
  #geom_smooth(method="lm", se=F) +
  facet_wrap(vars(strain), ncol=2)

# Conclusion: not a very strong link between traits and abundance, with some exceptions for the cyanos

## How well can we predict pcgr and trait change?-----
stats_result <- modelling(data=data, 
                          var_to_nest_by = c("strain", "treat"),
                          formulas=c("dT ~ density + trait", #+ density*mean
                                     "dT ~ trait",
                                     "pcgr ~ trait + density",# + density*mean
                                     "pcgr ~ density")) %>%
  select(-data)

#Now join these models to all the data and make the predictions
data_preds <- data %>%
  ungroup() %>%
  nest_by(species, strain, treat) %>%
  left_join(stats_result, by=c("strain", "treat"), multiple = "all") %>%
  mutate(predictions = list(cbind(data, pred=predict.lm(model, newdata=data)))) %>%
  rowwise() %>%
  mutate(AIC=AIC(model)) %>%  
  ungroup() %>%
  mutate(response = case_when(form=="dT ~ density + trait" ~ "trait change",
                              form=="dT ~ trait" ~ "trait change",
                              form=="pcgr ~ trait + density" ~ "growth",
                              form=="pcgr ~ density" ~ "growth")) %>%
  ungroup()
  
# Plot model fit and compare (AIC)
data_preds_synth <- data_preds %>%
  mutate(predictors = case_when(form=="dT ~ density + trait" ~ "both",
                                form=="dT ~ trait" ~ "focal",
                                form=="pcgr ~ trait + density" ~ "both",
                                form=="pcgr ~ density" ~ "focal")) %>%
  select(c("species", "strain", "treat", "response", "AIC", "predictors")) %>%
  pivot_wider(names_from = "predictors", values_from = "AIC") %>%
  mutate(delta_AIC = both - focal) #if <0 then both predict better

ggplot(data_preds_synth) + 
  scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=focal, y=both, col=as.factor(treat), pch=response) + 
  #aes(x=response, y=delta_AIC, col=as.factor(treat), pch=trait) + 
  #geom_jitter(width = 0.25) +
  geom_point() +
  #geom_hline(yintercept = 0, lty="dotted") +
  geom_abline(intercept = 0, slope=1) +
  facet_wrap(vars(strain), ncol=2, scales="free") + 
  labs(x = "AIC, trait or abundance only", y = "AIC, both predictors", 
       col="treatment", pch="response") + 
  geom_abline(intercept = -2, slope=1, lty="dotted") +
  geom_abline(intercept = 2, slope=1, lty="dotted")
  
#ggsave(paste0(model_system, "_", response, ".pdf"), 
#width=5, height = 4, device = "pdf")
ggsave(paste0(model_system, "_", response, ".pdf"), 
       width=5, height = 5, device = "pdf")

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

