#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT TOOLS AND DATA -----------------------
source("Tools and data.R")
# PICK DATA SOURCE ---------------
model_system <- "cyano" #cilia or cyano
data         <- get(paste("data_",model_system, sep=""))
# DO ANALYSES --------------------
## Just plot of dT vs. trait and pcgr vs. pop. ------
ggplot(data) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=density, y=pcgr, col=treat) +
  geom_point() +
  geom_smooth(method="lm", se=F, lwd=0.5) +
  facet_wrap(vars(strain), ncol=2, scales="free") 

ggsave(paste0("plots/", model_system,"_pcgr.pdf"), 
       width=5, height = 4, device = "pdf")

ggplot(data) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=trait, y=dT, col=treat) +
  geom_point() +
  geom_smooth(method="lm", se=F, lwd=0.5) +
  facet_wrap(vars(strain), ncol=2, scales="free") 

ggsave(paste0("plots/", model_system,"_dT.pdf"), 
       width=5, height = 4, device = "pdf")

###Result--------
#Ciliates: 
#negative dd, but different for different treatments; 
#negative td, not different for different treatments.
#Cyanos:
#dd (mostly negative); different for different treatments. 
#td (mostly negative); different for different treatments.

## Task 1: Predict next time point within condition: single vs. two variables-----
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
  aes(x=response, y=delta_AIC, col=as.factor(treat)) + 
  #aes(x=response, y=delta_AIC, col=as.factor(treat), pch=trait) + 
  geom_jitter(width = 0.25) +
  #geom_point() +
  geom_hline(yintercept = 2, lty="dotted") +
  geom_hline(yintercept = -2, lty="dotted") +
  geom_hline(yintercept = 0) +
  #geom_abline(intercept = 0, slope=1) +
  facet_wrap(vars(strain), ncol=2, scales="free") + 
  labs(y=expression(paste("AIC"[full],"-AIC"[single])), col="treatment")
  
#ggsave(paste0(model_system, "_", response, ".pdf"), 
#width=5, height = 4, device = "pdf")
ggsave(paste0("plots/", model_system,"_AIC.pdf"), 
       width=5, height = 4, device = "pdf")

# Plot model fits
data_preds <- data_preds %>%
  select(-c("data", "model")) %>%
  unnest(predictions)
#, form=="pcgr ~ density"
ggplot(data_preds %>% filter(response=="growth")) + 
  theme_bw() + 
  scale_shape_manual(values=0:10) + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=pcgr, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  geom_smooth(method="lm", se=F, lwd=0.5)+
  facet_grid(vars(form)) + #, scales="free"
  geom_abline(slope=1, intercept=0) 

ggsave(paste0("plots/", model_system,"growth.pdf"), 
       width=5, height = 4, device = "pdf")

ggplot(data_preds %>% filter(response=="trait change")) + 
  theme_bw() + 
  scale_shape_manual(values=0:10) + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=dT, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  geom_smooth(method="lm", se=F, lwd=0.5)+
  facet_grid(vars(form)) + #
  geom_abline(slope=1, intercept=0) 

ggsave(paste0("plots/", model_system,"trait.pdf"), 
       width=5, height = 4, device = "pdf")

###Result--------
#Ciliates: 
#Combining two predictors rarely helped to predict growth or trait change
#Cyanos:
#Combining two predictors rarely helped to predict growth, but did improve trait change prediction

## Task 2: Predict next time point across conditions: single vs. two predictors-----
stats_result <- modelling(data=data%>%filter(treat=="C"), 
                          var_to_nest_by = c("strain"),
                          formulas=c("dT ~ trait", 
                                     "dT ~ trait + density",
                                     "pcgr ~ density",
                                     "pcgr ~ trait + density")) %>%
  select(-data)

#Now join these models to all the data and make the predictions
data_preds <- data %>%
  ungroup() %>%
  nest_by(strain) %>%
  left_join(stats_result, by=c("strain"), multiple = "all") %>%
  mutate(response = case_when(length(grep("dT", form))>0 ~ "trait change",
                              length(grep("pcgr", form))>0 ~ "growth")) %>%
  mutate(predictor = case_when(length(grep("trait + density", form, fixed=T))>0 ~ "both",
                               length(grep("density", form))>0 ~ "one",
                               length(grep("trait", form))>0 ~ "one")) %>%
  #select(-form) %>%
  mutate(data_and_pred = list(cbind(data, prediction=predict.lm(object=model, newdata=data)))) %>%
  select(c("strain", "response", "form", "predictor", "data_and_pred")) %>%
  unnest(c("data_and_pred")) %>%
  #pivot_wider(values_from = "prediction", names_from = "response")
  mutate(error =  case_when(response=="trait change" ~ abs(prediction - dT),
                            response=="growth" ~ abs(prediction - pcgr)))

#Now plot, basic plot first: growth
ggplot(data_preds %>% filter(response=="growth")) + 
  aes(x=pcgr, y=prediction, col=treat, pch=strain) + #
  geom_point() +
  geom_smooth(method="lm", se=F, lwd=0.5) + 
  #scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  facet_grid(vars(form)) + #
  geom_abline(intercept=0, slope=1)

#Now plot trait change
ggplot(data_preds %>% filter(response=="trait change")) + 
  aes(x=dT, y=prediction, col=treat, pch=strain) + #
  geom_point() +
  geom_smooth(method="lm", se=F, lwd=0.5) + 
  #scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  facet_grid(vars(form)) + #
  geom_abline(intercept=0, slope=1)

#Now error plot
data_preds_synth <- data_preds %>%
  group_by(strain, response, predictor, treat) %>%
  summarise(error = sum(error)) %>%
  pivot_wider(names_from = predictor, values_from = error) %>%
  mutate(delta_error = (both - one)/one)

ggplot(data_preds_synth) + 
  scale_shape_manual(values=0:10) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=response, y=delta_error, col=as.factor(treat)) + 
  #aes(x=response, y=delta_AIC, col=as.factor(treat), pch=trait) + 
  geom_jitter(width = 0.25) +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(strain), ncol=2) + 
  labs(y=expression(paste("Error"[full],"-Error"[single])), 
       col="treatment")

###Result--------
#Ciliates: 
#Combining two predictors rarely helped to predict growth or trait change to other conditions
#Cyanos:
#Combining two predictors rarely helped to predict growth but could improve or worsen trait change prediction.

## Leftovers -------
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
  facet_wrap(vars(strain), ncol=2) +
  labs(col="treatment")

ggsave(paste0("plots/", model_system,"corr.pdf"), 
       width=5, height = 4, device = "pdf")

# Conclusion: not a very strong link between traits and abundance, with some exceptions for the cyanos

## Basic stats: treatment effects on dd of pcgr and td of dT -----
stats_result <- modelling(data=data, 
                          var_to_nest_by = c("strain", "treat"),
                          formulas=c("dT ~ trait")) %>%#"dT ~ trait", "pcgr ~ density"
  mutate(model_summary = 
           list(as_tibble(rownames_to_column(as.data.frame(summary(model)$coefficients),
                                             var = "predictor")))) %>%
  unnest(model_summary)  %>%
  rowwise() %>%
  mutate(type_of_pred = ifelse(length(grep("density", predictor)>0)|length(grep("trait", predictor)>0), "slope", "intercept")) %>%
  ungroup() %>%
  mutate(Estimate_sig = ifelse(`Pr(>|t|)`<0.05, Estimate, NA))

plot <- ggplot(stats_result) +#%>% filter(`Pr(>|t|)`<0.05)
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=strain, y=Estimate, col=treat) +
  geom_point(shape=1, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), 
                position=position_dodge(width=0.5), width=.2) +
  geom_point(aes(x=strain, y=Estimate_sig), 
             position=position_dodge(width=0.5), 
             pch="*", cex=5, show.legend = F) +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() + 
  facet_wrap(vars(type_of_pred), scales="free") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste("plots/td_general_",model_system,".pdf", sep=""), plot=plot, 
       width = 7, height = 7)
#ggsave(paste("plots/td_general_",model_system,".pdf", sep=""), plot=plot, 
#       width = 5, height = 3)