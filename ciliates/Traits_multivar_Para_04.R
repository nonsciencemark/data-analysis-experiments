
library(tidyverse)
library(ggplot2)
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

Traj.avg <- read.csv("OverviewTraitsFull.csv", header = T) %>%
  filter(ID_spec == "Para_4")


Expl_vars <- Traj.avg[c("Temp", "Atrazine", "Days_fromstart", "ID_spec")]
Obs_vars <- Traj.avg[c(8:16)]

Traj.avg

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
pairs(Obs_vars[c(1:9)], lower.panel = panel.cor) # Keep all for now

Obs_vars.select <- Obs_vars

pairs(Obs_vars.select, lower.panel = panel.cor)

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
                         #"duration",
                         "mean_linearity",
                         "Parts_permL"
                         )]

Obs_sd <- Obs_vars.select[c("sd_ar", 
                            "sd_area", 
                            #"mean_ar",
                            "sd_speed",
                            #"duration",
                            "sd_linearity"
)]
Expl <- Expl_vars

#First column of the dataset as row labels: NOT USEFULL HERE AS STATION DONT RESULT IN UNIQUE ROWNAMES

Total <- cbind(Obs, Expl)
Totalsd <- cbind(Obs_sd, Total) %>%
  mutate(linearity_uppersd = mean_linearity + sd_linearity,
         linearity_lowersd = mean_linearity - sd_linearity,
         area_uppersd = mean_area + sd_area,
         area_lowersd = mean_area - sd_area,
         ar_uppersd = mean_ar + sd_ar,
         ar_lowersd = mean_ar - sd_ar,
         speed_uppersd = mean_speed + sd_speed,
         speed_lowersd = mean_speed - sd_speed) #%>%
  #select(-c(sd_linearity, sd_area, sd_ar, sd_speed))



#Do a PCA on mean traits and abundance
res.pca <- PCA(Total[c("mean_ar", 
                       "mean_area", 
                       #"mean_ar",
                       "mean_speed",
                       #"duration",
                       "mean_linearity",
                       "Parts_permL"
                       )])

summary(res.pca)

#Do a PCA on trait variance
res.pca2 <- PCA(Totalsd[c("sd_ar", 
                        "sd_area", 
                        #"mean_ar",
                        "sd_speed",
                        #"duration",
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


corrplot(varsd$cos2,method = "number",cl.lim = c(0,1),tl.col = "black",col=colorRampPalette(c("white","snow3","black"))(200))


#The sum of the representation of our variables on the different dimensions (here 2)
fviz_cos2(res.pca, choice = "var", axes = 1:2,ylim=c(0,1),fill = "lightgray", color = "black")+
  geom_text(aes(label=round(cos2, digits = 2)),fontface=2, vjust=1.6, color="black", size=3.5)+
  ggtitle(label = "Representation of our variables on the PCA dimensions 1 to 2")+
  theme(plot.title = element_text(hjust = 0.5))        ## 

fviz_cos2(res.pca2, choice = "var", axes = 1:2,ylim=c(0,1),fill = "lightgray", color = "black")+
  geom_text(aes(label=round(cos2, digits = 2)),fontface=2, vjust=1.6, color="black", size=3.5)+
  ggtitle(label = "Representation of our variables on the PCA dimensions 1 to 2")+
  theme(plot.title = element_text(hjust = 0.5)) 

#COrrelation circles: speed is the most important, parts per mL and aspect ratio contribute the least
fviz_pca_var(res.pca,axes = 1:2, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_blank() ) # logPO4, Chla, MLD/DCM are representing reality well, logNP and logNH4 badly

#fviz_pca_var(res.pca,axes = 2:3, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
#)+theme(plot.title = element_blank() ) # logNP, logN2_zavgsqrtPar and MLD/DCM represent reality well, rest badly


# variance in linearity is most important, variance in aspect ratio contributes the least
fviz_pca_var(res.pca2,axes = 1:2, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
)+theme(plot.title = element_blank() ) # logPO4, Chla, MLD/DCM are representing reality well, logNP and logNH4 badly

#fviz_pca_var(res.pca2,axes = 2:3, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE # Avoid text overlapping
#)+theme(plot.title = element_blank() ) # logNP, logN2_zavgsqrtPar and MLD/DCM represent reality well, rest badly

# axes 1 + 2 Temp: slight shift towards slower, smaller and non-linear moving cells
fviz_pca_biplot(res.pca, pointsize = 5, 
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


## axes 1,2, Culture age: more variation in speed over time
fviz_pca_biplot(res.pca2,
                pointsize = 5,# Variables color
                label="var",
                col.var = "black",
                pointshape = 21, fill.ind = Totalsd$Days_fromstart,
                legend.title = "Days from start",
                gradient.cols = c("yellow", "purple")
)+theme(plot.title = element_text(hjust = 0.5))


# Intermediate concludion: Culture age has a stronger effect than temp or atra

# Try an RDA with Temp, Atrazine and Days_fromstart:
library(vegan)

# normalize data:
Total <- Total %>%
  mutate(Atrazine = as.numeric(as.character(Atrazine)),
         Temp = as.numeric(as.character(Temp)))

Totalsd <- Totalsd %>%
  mutate(Atrazine = as.numeric(as.character(Atrazine)),
         Temp = as.numeric(as.character(Temp)))

m<-apply(Total[-c(9)],2,mean) #Calc mean for all columns
m
s<-apply(Total[-c(9)],2,sd) # idem for sd
s

Obsnorm<-as.data.frame(scale(Total[-c(9)],m,s)) #normalized data
Obsnorm$ID_spec <- Total$ID_spec
Obsnorm$Temp <- Total$Temp
Obsnorm$Atrazine <- Total$Atrazine
Obsnorm$Days_fromstart <- Total$Days_fromstart

m<-apply(Totalsd[-c(8)],2,mean) #Calc mean for all columns
m
s<-apply(Totalsd[-c(8)],2,sd) # idem for sd
s

Obsnormsd<-as.data.frame(scale(Totalsd[-c(8)],m,s)) #normalized data
Obsnormsd$ID_spec <- Totalsd$ID_spec
Obsnormsd$Temp <- Totalsd$Temp
Obsnormsd$Atrazine <- Totalsd$Atrazine
Obsnormsd$Days_fromstart <- Totalsd$Days_fromstart
#normalize expl data 

# put together:

TotalNorm <- cbind(Obsnorm, Obsnormsd)


# plot avg of traits with upper and lower sd values:
ggplot(data = Totalsd, mapping = aes(x = Days_fromstart, y = mean_linearity, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_ribbon(data = Totalsd, mapping = aes(x = Days_fromstart, ymin = mean_linearity - sd_linearity, ymax = mean_linearity + sd_linearity,
                                            fill = as_factor(Temp)), alpha = 1/10) +
  geom_smooth(se=F) +
  labs(title = "Mean Movement linearity over time") +
  theme_classic() 

ggplot(data = Totalsd, mapping = aes(x = Days_fromstart, y = mean_area, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_ribbon(data = Totalsd, mapping = aes(x = Days_fromstart, ymin = mean_area - sd_area, ymax = mean_area + sd_area,
                                            fill = as_factor(Temp)), alpha = 1/10) +
  geom_smooth(se=F) +
  labs(title = "Cell size over time") +
  theme_classic() 

ggplot(data = Totalsd, mapping = aes(x = Days_fromstart, y = mean_ar, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_ribbon(data = Totalsd, mapping = aes(x = Days_fromstart, ymin = mean_ar - sd_ar, ymax = mean_ar + sd_ar,
                                            fill = as_factor(Temp)), alpha = 1/10) +
  geom_smooth(se=F) +
  labs(title = "Cell shape over time") +
  theme_classic() 

ggplot(data = Totalsd, mapping = aes(x = Days_fromstart, y = mean_speed, color = as_factor(Temp))) +
  facet_wrap(facets = ~Totalsd$Atrazine) +
  geom_point() +
  geom_ribbon(data = Totalsd, mapping = aes(x = Days_fromstart, ymin = mean_speed - sd_speed, 
                                            ymax = mean_speed + sd_speed,
                                            fill = as_factor(Temp)), alpha = 1/10) +
  geom_smooth(se=F) +
  labs(title = "Movement speed over time") +
  theme_classic() 

