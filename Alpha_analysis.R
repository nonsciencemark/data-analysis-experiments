library(tidyverse)
library(lubridate)
library(mgcv)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")#
#devtools::install_github("ctkremer/growthTools")
#library(growthTools)

# IMPORT AND MAKE UNIFORM THE DATA -----------------------
source("Streamline data.R")
# PICK DATA SOURCE ---------------
model_system <- "cyano" #cilia or cyano
data         <- get(paste("data_",model_system, sep=""))
# DO ANALYSES --------------------
##check if dd changes with treatment -------
ggplot(data) +
  scale_colour_manual(values=c("black", cbPalette)) + 
  aes(x=density,y=pcgr, col=as.factor(atrazine)) +
  geom_point() +
  facet_wrap(vars(strain, temperature), scales="free", 
             ncol=length(unique(data$temperature)),#,
             labeller = label_bquote(paste("T=", .(temperature),
                                           ", strain=", .(strain)))) +
  geom_smooth(method=lm, aes(x=density, y=pcgr, col=as.factor(atrazine)),
              formula=y ~ poly(x,1), se=F) + 
  labs(x="Inds per mL", y="pcgr", col="atrazine")
ggsave(paste("dd_",model_system,".pdf", sep=""), width=1+2*length(unique(data$temperature)), 
       height = 4*length(unique(data$temperature)), device = "pdf")
model <- lm(pcgr~(Parts_permL+as.factor(Atrazine)+as.factor(Temp)+strain)^2, 
            data=data_cilia)

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
  ggsave(paste("dd_", model_system, "_",trait_i,".pdf", sep=""), plot=plot_i, 
         width = 1+2*length(unique(data$temperature)), 
         height = 4*length(unique(data$temperature)))
  model_i <- lm(log10(mean)~(log10(density)+as.factor(atrazine)+as.factor(temperature)+strain)^2, 
              data=subset(data, trait==trait_i))
  print(trait_i)
  print(summary(model_i))
  }

## fit a reference model of dd and trait dependence on density--------
# This model only uses the control data
data_ref <- data %>% 
  filter(atrazine%in%c("0","no"), temperature %in%c("normal", "20")) 
###first dd----------------
ref_model_dd <- lm(data_ref, formula = pcgr ~ poly(density,1)*strain + strain)

###then traits----------------
ref_model_trait <- lm(data_ref,  formula = mean ~ poly(log10(density),1)*strain*trait + 
                        strain*trait + strain + trait)
###add predictions to the original data frame------------
data$pcgr_ref <- predict.lm(ref_model_dd, newdata = data)
data$mean_ref <- predict.lm(ref_model_trait, newdata = data)
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
  ggsave(paste("delta_",model_system,"_",trait_i,".pdf", sep=""), plot=plot_i, 
         width = 1+2*length(unique(data$temperature)), 
         height = 4*length(unique(data$temperature)))
  model_i <- lm(delta_trait~(delta_pcgr+as.factor(atrazine)+as.factor(temperature)+strain)^2, 
                data=subset(data, trait==trait_i))
  print(trait_i)
  print(summary(model_i))
}


