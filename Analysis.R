#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT TOOLS AND DATA -----------------------
source("Tools and data.R")
# PICK DATA SOURCE ---------------
model_system <- "cyano" #cilia or cyano
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
  nest_by(strain, trait, treat) %>%
  left_join(stats_result, by=c("strain", "trait", "treat"), multiple = "all") %>%
  mutate(predictions = list(cbind(data, pred=predict.lm(model, newdata=data)))) %>%
  select(-c("data", "model")) %>%
  unnest(predictions) %>%
  rowwise() %>%
  mutate(response = list(ifelse(length(grep("dT", form))==0, 
                                "growth", "trait")), .after="form")

ggplot(data_preds %>% filter(response=="growth")) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=pcgr, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  facet_grid(form~trait) + #, scales="free"
  geom_abline(slope=1, intercept=0) 

ggplot(data_preds %>% filter(response=="trait")) + 
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=dT, y=pred, col=as.factor(treat), pch=strain) + 
  geom_point() + 
  #geom_smooth(method="lm", se=F)+
  facet_wrap(trait~form, scales="free") + #
  geom_abline(slope=1, intercept=0) 
