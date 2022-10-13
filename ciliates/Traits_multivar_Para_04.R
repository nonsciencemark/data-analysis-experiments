
library(tidyverse)
library(lubridate)
library(mgcv)
'''
Traj.avg <- read.csv("U:/DIVERCE_experiments/DIVERCE_WP01_exp01/Analysis/Overview_dataframe/OverviewTraitsFull.csv", header = T) %>%
  mutate(Parts_permL = Count/830*1000)

for (i in 1:length(Traj.avg$ID_spec)) {
  if(Traj.avg$ID_spec[i] == "Tetra_1") {
    Traj.avg$Parts_permL[i] <- Traj.avg$Parts_cor[i]
    } else {
      if(Traj.avg$ID_spec[i] == "Tetra_2") {
        Traj.avg$Parts_permL[i] <- Traj.avg$Parts_cor[i]
      } else {
        Traj.avg$Parts_permL[i] <- Traj.avg$Parts_permL[i]
      }
      
    }
}

write.csv(Traj.avg, file = 'U:/DIVERCE_experiments/DIVERCE_WP01_exp01/Analysis/Overview_dataframe/OverviewTraitsFull.csv', 
          row.names = F)
'''

# CILIATES ------
## Import data and compute pcgr, and only keep specific species-trait combinations----
data_cilia <- read.csv("ciliates/OverviewTraitsFull.csv") %>%
  #filter(ID_spec == "Para_4") %>%
  #rename(strain = ID_spec) %>%
  group_by(Atrazine, Temp, ID_spec) %>%
  mutate(pcgr = log(lead(Parts_permL, 1)/Parts_permL)/(lead(Days_fromstart, 1)-Days_fromstart),
         .after = Parts_permL) %>% 
  filter(pcgr > -5) %>% #Kick out extremely negative pcgr
  ungroup() %>%
  mutate(phase = cut_interval(log10(Parts_permL), n=5, labels = FALSE)) %>%
  select(ID_spec, Sample,  Days_fromstart, Treatment, Temp, Atrazine, contains(c("mean_", "sd")), pcgr, Parts_permL)
  #pivot_longer(mean_ar:sd_linearity, names_to="trait") %>%
  #separate(trait, into=c("stat", "trait"), sep="_") %>%
  #pivot_wider(names_from=stat, values_from=value) %>%
  #mutate(cv=sd/mean) %>%
  
  #filter(((species=="Spiro")&(trait=="linearity"))|((species=="Tetra")&(trait=="area"))
   #      |((species=="Para")&(trait=="speed"))|((species=="Loxo")&(trait=="linearity")))


Expl_vars <- data_cilia[c("Temp", "Atrazine", "Days_fromstart", "ID_spec")]
Obs_vars <- data_cilia[c(7:16)]

### FROM KRISTINA MOJICA:
#Function for multi-panel Cleveland dotplot.
#The input file must contain no categorical variables
Mydotplot <- function(DataSelected){
  
  P <- dotplot(as.matrix(as.matrix(DataSelected)),
               groups=FALSE,
               strip = strip.custom(bg = 'white',
                                    par.strip.text = list(cex = 1.2)),
               scales = list(x = list(relation = "free", draw = TRUE),
                             y = list(relation = "free", draw = FALSE)),
               col=1, cex  = 0.5, pch = 16,
               xlab = list(label = "Value of the variable", cex = 1.5),
               ylab = list(label = "Order of the data from text file", cex = 1.5))
  
  print(P)  
}

library(car)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#pairs correlation plot
pairs(Obs_vars[c(1:10)], lower.panel = panel.cor) # Keep all for now

Obs_vars.select <- Obs_vars

library(Hmisc)
hist.data.frame(Obs_vars.select) # not much use actually as data hasnt been standardized yet..
boxplot(Obs_vars.select)



###########NOW WITH NILS' WAY OF PERFORMING PCA#########################################

#My code
library(FactoMineR)
library(factoextra)


#Make same selection as in normalized2 but denormalized cause pca in facominer already does that:
Obs <- Obs_vars.select[c("mean_ar", 
                         "mean_area", 
                         #"mean_ar",
                         "mean_speed",
                         "pcgr",
                         "mean_linearity",
                         "Parts_permL"
                         )]

Obs_sd <- Obs_vars.select[c("sd_ar", 
                            "sd_area", 
                            "pcgr",
                            "sd_speed",
                            #"duration",
                            "sd_linearity"
)]
Expl <- Expl_vars

#First column of the dataset as row labels: NOT USEFULL HERE AS STATION DONT RESULT IN UNIQUE ROWNAMES

Total <- cbind(Obs, Expl)
Totalsd <- full_join(Obs_sd, Total)



#Do a PCA on mean traits and abundance
res.pca <- PCA(Total[c("mean_ar", 
                       "mean_area", 
                       #"mean_ar",
                       "mean_speed",
                       "mean_linearity"
                       )])

summary(res.pca)

#Do a PCA on trait variance
res.pca2 <- PCA(Totalsd[c("sd_ar", 
                          "sd_area", 
                          "sd_speed",
                          "sd_linearity"
)])

summary(res.pca2)


var <- get_pca_var(res.pca)

varsd <- get_pca_var(res.pca2)

#The weight of the different dimensions of the PCA
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70))+
  ggtitle(label = "Representation of our dataset on the PCA dimensions")+
  theme(plot.title = element_text(hjust = 0.5))

fviz_eig(res.pca2, addlabels = TRUE, ylim = c(0, 70))+
  ggtitle(label = "Representation of our dataset on the PCA dimensions")+
  theme(plot.title = element_text(hjust = 0.5))

#--> First 3? axes seem relevant

#The weight of the different variables on the different dimensions
library("corrplot")
corrplot(var$cos2,method = "number",cl.lim = c(0,1),tl.col = "black",col=colorRampPalette(c("white","snow3","black"))(200))
# pcgr is represented in dim 1 and 3

corrplot(varsd$cos2,method = "number",cl.lim = c(0,1),tl.col = "black",col=colorRampPalette(c("white","snow3","black"))(200))


#The sum of the representation of our variables on the different dimensions (here 2)
fviz_cos2(res.pca, choice = "var", axes = 1:3,ylim=c(0,1),fill = "lightgray", color = "black")+
  geom_text(aes(label=round(cos2, digits = 2)),fontface=2, vjust=1.6, color="black", size=3.5)+
  ggtitle(label = "Representation of our variables on the PCA dimensions 1 to 3")+
  theme(plot.title = element_text(hjust = 0.5))        ## 

fviz_cos2(res.pca2, choice = "var", axes = 1:3,ylim=c(0,1),fill = "lightgray", color = "black")+
  geom_text(aes(label=round(cos2, digits = 2)),fontface=2, vjust=1.6, color="black", size=3.5)+
  ggtitle(label = "Representation of our variables on the PCA dimensions 1 to 3")+
  theme(plot.title = element_text(hjust = 0.5)) 

#COrrelation circles: speed is the most important, parts per mL and aspect ratio contribute the least
fviz_pca_var(res.pca,axes = 1:2, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_blank() ) # logPO4, Chla, MLD/DCM are representing reality well, logNP and logNH4 badly

fviz_pca_var(res.pca,axes = c(1,3), col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_text()) + # logNP, logN2_zavgsqrtPar and MLD/DCM represent reality well, rest badly
labs(title = "Axes 1 vs 3")

# variance in linearity is most important, variance in aspect ratio contributes the least
fviz_pca_var(res.pca2,axes = 1:2, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_blank() ) # logPO4, Chla, MLD/DCM are representing reality well, logNP and logNH4 badly

fviz_pca_var(res.pca2,axes = 2:3, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_blank() ) # logNP, logN2_zavgsqrtPar and MLD/DCM represent reality well, rest badly

# axes 1 + 2 Temp: slight shift towards slower, smaller and non-linear moving cells
fviz_pca_biplot(res.pca, pointsize = 5, 
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Total$ID_spec),
                legend.title = "Temperature",
                addEllipses = TRUE,
                repel = TRUE)

fviz_pca_biplot(res.pca, pointsize = 5, 
                axes = c(1,3),
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Total$Temp),
                legend.title = "Temperature",
                addEllipses = TRUE,
                repel = TRUE)

# less variance with increasing temp? meaning narrower niche?
fviz_pca_biplot(res.pca2, pointsize = 5, 
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Totalsd$ID_spec),
                legend.title = "Temperature",
                addEllipses = TRUE,
                repel = TRUE)

fviz_pca_biplot(res.pca2, pointsize = 5,
                axes = c(1,3),
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Totalsd$Temp),
                legend.title = "Temperature",
                addEllipses = TRUE,
                repel = TRUE)


# axes 1 + 2 Atra: No effect that I see
fviz_pca_biplot(res.pca, pointsize = 5, 
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Total$Atrazine),
                legend.title = "Atrazine",
                addEllipses = TRUE,
                repel = TRUE)

fviz_pca_biplot(res.pca, pointsize = 5, 
                axes = c(1,3),
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Total$Atrazine),
                legend.title = "Atrazine",
                addEllipses = TRUE,
                repel = TRUE)

# axes 1 + 2 Atra: tiny bit less variance with increasing temp?
fviz_pca_biplot(res.pca2, pointsize = 5, 
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = as.factor(Totalsd$Atrazine),
                legend.title = "Atrazine",
                addEllipses = TRUE,
                repel = TRUE)

## axes 1,2, Culture age: smaller cells over time
fviz_pca_biplot(res.pca,
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Total$Days_fromstart,
                legend.title = "Days from start",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))

fviz_pca_biplot(res.pca,
                axes = c(1,3),
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Total$Days_fromstart,
                legend.title = "Days from start",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))

## axes 1,2, Culture age: more variation in speed over time
fviz_pca_biplot(res.pca2,
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Totalsd$Days_fromstart,
                legend.title = "Days from start",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))


# pcgr as color gradient:

## axes 1,2, 
fviz_pca_biplot(res.pca,
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Total$pcgr,
                legend.title = "Per Capita Growth rate",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))

fviz_pca_biplot(res.pca,
                axes = c(1,3),
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Total$pcgr,
                legend.title = "Per Capita Growth rate",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))

## axes 1,2, 
fviz_pca_biplot(res.pca2,
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Totalsd$pcgr,
                legend.title = "Per Capita Growth rate",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))



# plot pcgr vs traits to see for correlations:
ggplot(data = Totalsd, mapping = aes(x = pcgr, y = mean_linearity, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine + ID_spec) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "Mean Movement linearity versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = mean_area, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "Cell size versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = mean_ar, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "Cell shape versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = mean_speed, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "Movement speed versus growth") +
  theme_classic(base_size = 24) 


# plot avg of traits with upper and lower sd values:
ggplot(data = Totalsd, mapping = aes(x = pcgr, y = sd_linearity, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "sd Movement linearity versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = sd_area, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "sd Cell size versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = sd_ar, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "sd Cell shape versus growth") +
  theme_classic(base_size = 24) 

ggplot(data = Totalsd, mapping = aes(x = pcgr, y = sd_speed, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_smooth(se=T, method = "lm") +
  labs(title = "sd Movement speed versus growth") +
  theme_classic(base_size = 24) 


# part 1.5 # CODE FROM QZ
# mean trait correlation,Para_4
unique(Totalsd$ID_spec)
data2 <-Totalsd

cor(data2$mean_ar, data2$mean_area, method = c("pearson"))  #
cor(data2$mean_ar, data2$mean_speed, method = c("pearson")) #
cor(data2$mean_ar, data2$mean_linearity, method = c("pearson")) #

cor(data2$mean_area, data2$mean_speed, method = c("pearson")) #
cor(data2$mean_area, data2$mean_linearity, method = c("pearson")) #

cor(data2$mean_speed, data2$mean_linearity, method = c("pearson")) #0.5464877
a2<-summary(lm(data2$mean_linearity~data2$mean_speed))
a2

data2$Atrazine <- as.factor(data2$Atrazine)


library("ggpubr")
library("ggplot2")
a<-ggscatter(data2, x = "pcgr", y = "mean_speed",
             shape = "ID_spec",
             color = "ID_spec", 
             size=3,
             add = "reg.line", conf.int = FALSE,
             cor.coef = TRUE, cor.method = "pearson",
             title= "",  
             xlab = "pcgr", ylab = "mean_speed"
)

a+theme_bw()+
  facet_wrap(~Atrazine + Temp) +
  #geom_abline(intercept =a2$coefficients[1,1], slope = a2$coefficients[2,1], color="black", linetype="solid", size=1)+
  theme(plot.title = element_text(hjust = 0.5,size=15), legend.position='right')
#geom_ab

a<-ggscatter(data2, x = "pcgr", y = "mean_linearity",
             shape = "ID_spec",
             color = "ID_spec", 
             size=3,
             add = "reg.line", conf.int = FALSE,
             cor.coef = TRUE, cor.method = "pearson",
             title= "",  
             xlab = "pcgr", ylab = "mean_linearity"
)

a+theme_bw()+
  facet_wrap(~Atrazine + Temp) +
  #geom_abline(intercept =a2$coefficients[1,1], slope = a2$coefficients[2,1], color="black", linetype="solid", size=1)+
  theme(plot.title = element_text(hjust = 0.5,size=15), legend.position='right')
#geom_ab

a<-ggscatter(data2, x = "pcgr", y = "mean_area",
             shape = "ID_spec",
             color = "ID_spec", 
             size=3,
             add = "reg.line", conf.int = FALSE,
             cor.coef = TRUE, cor.method = "pearson",
             title= "",  
             xlab = "pcgr", ylab = "mean_area"
)

a+theme_bw()+
  facet_wrap(~Atrazine + Temp) +
  #geom_abline(intercept =a2$coefficients[1,1], slope = a2$coefficients[2,1], color="black", linetype="solid", size=1)+
  theme(plot.title = element_text(hjust = 0.5,size=15), legend.position='right')
#geom_ab

a<-ggscatter(data2, x = "pcgr", y = "mean_ar",
             shape = "ID_spec",
             color = "ID_spec", 
             size=3,
             add = "reg.line", conf.int = FALSE,
             cor.coef = TRUE, cor.method = "pearson",
             title= "",  
             xlab = "pcgr", ylab = "mean_aspect_ratio"
)

a+theme_bw()+
  facet_wrap(~Atrazine + Temp) +
  #geom_abline(intercept =a2$coefficients[1,1], slope = a2$coefficients[2,1], color="black", linetype="solid", size=1)+
  theme(plot.title = element_text(hjust = 0.5,size=15), legend.position='right')
#geom_ab