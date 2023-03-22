library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
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
  filter(atrazine%in%c("0","no"), temperature %in%c("normal", "20")) 
stats_result_r <- modelling(data=data_ref, 
                          var_to_nest_by = c("strain", "trait"),
                          formula="pcgr ~ density + mean + density*mean")
stats_result_t <- modelling(data=data_ref, 
                            var_to_nest_by = c("strain", "trait"),
                            formula="dT ~ density + mean + density*mean")
#Now join these models to all the data (just the coefficients suffice)

n <- length(unique(data$strain))
plot_dd <- ggplot(stats_result %>% filter(`Pr(>|t|)`<0.05/n)) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  aes(x=predictor, y=Estimate, col=strain) +
  geom_point(shape=1, size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), 
                position=position_dodge(width=0.5), width=.2) +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() + 
  facet_wrap(vars(type_of_pred), scales="free", ncol=1)

ggsave(paste("plots/dd_general_",model_system,".pdf", sep=""), plot=plot_dd, 
       width = 5, height = 3)

##check if dependence of trait on density changes with treatment -------
for (trait_i in unique(data$trait)){
  plot_i <- ggplot(subset(data, trait==trait_i)) +
  scale_colour_manual(values=c("black", cbPalette)) + 
  aes(x=log10(density),y=log10(mean), col=as.factor(atrazine)) +
  geom_point() +
  facet_wrap(vars(strain, temperature), scales="free", 
             ncol=length(unique(data$temperature)),#,
             labeller = label_bquote(paste("T=", .(temperature),
                                           ", strain=", .(strain)))) +
  geom_smooth(method=lm, aes(x=log10(density), y=log10(mean), 
                             col=as.factor(atrazine)),
              formula=y ~ poly(x,1), se=F)+
  labs(x="log10(density)", y="log10(trait)", col="atrazine")
  ggsave(paste("plots/dd_", model_system, "_",trait_i,".pdf", sep=""), plot=plot_i, 
         width = 1+2*length(unique(data$temperature)), 
         height = 4*length(unique(data$temperature)))
}

### do the stats --------
stats_result <- modelling(data=data, 
                          var_to_nest_by = c("strain", "trait"),
                          formula="log10(mean)~ atrazine + temperature + log10(density)*atrazine + log10(density)*temperature")%>%
  rowwise() %>%
  mutate(type_of_pred = ifelse(length(grep("density", predictor)>0), "slope", "intercept")) %>%
  ungroup()
n <- length(unique(data$strain)) * length(unique(data$trait))

plot_dd <- ggplot(stats_result %>% filter(`Pr(>|t|)`<0.05/n)) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  scale_shape_manual(values=c(0:10)) +
  aes(x=predictor, y=Estimate, col=strain, shape=trait) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), 
                width=.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() + 
  facet_wrap(vars(type_of_pred), scales="free", ncol=1)

ggsave(paste("plots/dd_general_trait_",model_system,".pdf", sep=""), plot=plot_dd, 
       width = 6, height = 4)

## fit a reference model of dd and trait dependence on density--------
# This model only uses the control data
###first dd----------------
ref_model_dd <- lm(data_ref, formula = pcgr ~ poly(density,1)*strain + strain)

###then traits----------------
ref_model_trait <- lm(data_ref,  formula = mean ~ poly(log10(density),1)*strain*trait + 
                        strain*trait + strain + trait)
###add predictions to the original data frame------------
data$pcgr_ref <- predict.lm(ref_model_dd, newdata = data)#predictions of pcgr
data$mean_ref <- predict.lm(ref_model_trait, newdata = data)#predictions of mean trait
data <- data %>%
  mutate(delta_pcgr = pcgr-pcgr_ref) %>%
  mutate(delta_trait = mean-mean_ref)

##now check if delta with ref model depends on the effect on trait----------
for (trait_i in unique(data$trait)){
  plot_i <- ggplot(subset(data, trait==trait_i)) +
    theme_bw() + 
    scale_colour_manual(values=c("black", cbPalette)) + 
    geom_point(aes(x=delta_trait, y=delta_pcgr, col=as.factor(atrazine)), 
               alpha=0.8) + 
    geom_smooth(method=lm, aes(x=delta_trait, y=delta_pcgr, 
                               col=as.factor(atrazine)),
                formula=y ~ poly(x,1), se=F) + #,  col="black"
    facet_wrap(vars(strain, temperature), scales="free", 
               ncol=length(unique(data$temperature)),#,
               labeller = label_bquote(paste("T=", .(temperature),
                                             ", strain=", .(strain)))) +
    labs(x=expression(paste(delta[trait])), y=expression(paste(delta[pcgr])), col="atrazine")
  ggsave(paste("plots/delta_",model_system,"_",trait_i,".pdf", sep=""), plot=plot_i, 
         width = 1+2*length(unique(data$temperature)), 
         height = 4*length(unique(data$temperature)))
}

### do the stats --------
stats_result <- modelling(data=data, 
                          var_to_nest_by = c("strain", "trait"),
                          formula="delta_pcgr~ atrazine + temperature + delta_trait*atrazine + delta_trait*temperature")%>%
  rowwise() %>%
  mutate(type_of_pred = ifelse(length(grep("delta", predictor)>0), "slope", "intercept")) %>%
  ungroup()
n <- length(unique(data$strain)) * length(unique(data$trait))

plot_vs <- ggplot(stats_result %>% filter(`Pr(>|t|)`<0.05/n)) +
  theme_bw() + 
  scale_colour_manual(values=cbPalette) + 
  scale_shape_manual(values=c(0:10)) +
  aes(x=predictor, y=Estimate, col=strain, shape=trait) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), 
                width=.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() + 
  facet_wrap(vars(type_of_pred), scales="free", ncol=1)

ggsave(paste("plots/delta_vs_delta_",model_system,".pdf", sep=""), plot=plot_vs, 
       width = 5, height = 4)


