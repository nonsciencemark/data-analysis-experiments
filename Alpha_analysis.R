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
  aes(x=density,y=pcgr, col=atrazine) +
  geom_point() +
  facet_wrap(vars(strain, temperature), scales="free", 
             ncol=length(unique(data$temperature)),#,
             labeller = label_bquote(paste("T=", .(as.character(temperature)),
                                           ", strain=", .(as.character(strain))))) +
  geom_smooth(method=lm, aes(x=density, y=pcgr, col=as.factor(atrazine)),
              formula=y ~ poly(x,1), se=F) + 
  labs(x="Inds per mL", y="pcgr", col="atrazine")
ggsave(paste("dd_",model_system,".pdf", sep=""), width=1+2*length(unique(data$temperature)), 
       height = 4*length(unique(data$temperature)), device = "pdf")
### do the stats --------
n <- length(unique(data$strain))
modelling <- data %>%
  ungroup()%>%
  nest_by(strain) %>%
  mutate(model = list(summary(lm(pcgr~ atrazine + temperature + density*atrazine + density*temperature, 
                                 data = data))$coefficients)) %>%
  mutate(test = list((model[which(model[,"Pr(>|t|)"]<0.05/n),])))

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
}

### do the stats --------
n <- length(unique(data$strain)) * length(unique(data$trait))
modelling <- data %>%
  ungroup()%>%
  nest_by(strain, trait) %>% 
  mutate(model = list(summary(lm(log10(mean)~ atrazine + temperature + log10(density)*atrazine + log10(density)*temperature, 
                                 data = data))$coefficients)) %>%
  mutate(test = list((as.data.frame(model)[which(model[,"Pr(>|t|)"]<0.05/n),])))

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
  ggsave(paste("delta_",model_system,"_",trait_i,".pdf", sep=""), plot=plot_i, 
         width = 1+2*length(unique(data$temperature)), 
         height = 4*length(unique(data$temperature)))
}

### do the stats --------
n <- length(unique(data$strain)) * length(unique(data$trait))
modelling <- data %>%
  ungroup()%>%
  nest_by(strain, trait) %>%
  mutate(model = list(summary(lm(delta_pcgr ~ atrazine + temperature + 
                                   delta_trait*atrazine + delta_trait*temperature, 
                                 data = data))$coefficients)) %>%
  mutate(test = list((as.data.frame(model)[which(model[,"Pr(>|t|)"]<0.05/n),])))




