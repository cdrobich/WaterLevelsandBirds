
Species <- read.csv("Data/Species matrix_column relativized.csv")
sp.env <- read.csv("Data/Env_variables.csv")


library(vegan)
library(ggplot2)
library(gridExtra)
library(plyr)
library(ggrepel)
library(tidyverse)

Species <- Species[,2:24]

set.seed(105)  #reproduce results (used the time)

species.nmds <- metaMDS(Species, distance = "bray", k = 3, autotransform = FALSE, trymax = 100)
species.nmds

stressplot(species.nmds)
species.nmds$iters #124 iterations

#Data:     Species 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1555859 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘Species’ 

ordiplot(species.nmds, choices = c(1,2), type = "points") # axis 1, 2
ordiplot(species.nmds, choices = c(1,3), type = "points") # axis 1, 3


# extract the scores into a data frame

data.scores <- as.data.frame(scores(species.nmds)) 

data.scores$site <- sp.env$Site

data.scores$Water.depth <- sp.env$Depth
data.scores$Year <- as.factor(sp.env$Year) 
data.scores$vegetation <- sp.env$VegType
data.scores$vegYr <- sp.env$VegYr

str(data.scores)

#looking at what taxa are correlated w/ sites
alltaxa <- envfit(species.nmds, Species, permutations = 999, choices = c(1,2,3)) 
alltaxa


speciesfit12 <- Species[1]
speciesfit12 #my lazy way of adding which taxa I want to appear on the plot, usually r2 > .20
speciesfit12$AMBI <- Species$AMBI  ## r2 = 0.2118
speciesfit12$MAWR <- Species$MAWR  ## r2 = 0.5221
speciesfit12$SORA <- Species$SORA  ## r2 = 0.2732
speciesfit12$SWSP <- Species$SWSP  ##r2 = 0.3056
speciesfit12$VIRA <- Species$VIRA ## r2 = 0.2639
speciesfit12$AMWO <- Species$AMWO  ## r2 = 0.2171
speciesfit12$RWBL <- Species$RWBL  ##r2 = 0.3068
speciesfit12$YWAR <- Species$YWAR ##r2 = 0.3790
speciesfit12$SOSP <- Species$SOSP  ##r2 = 0.5300
speciesfit12$NOCA <- Species$NOCA  ## r2 = 0.2018
speciesfit12$TRES <- Species$TRES  ##r2 = 0.5408
speciesfit12$BARS <- Species$BARS  ##r2 = 0.3619
speciesfit12$BANS <- Species$BANS ## r2 = 0.2126
speciesfit12$EAKI <- Species$EAKI ##r2 = 0.2882
speciesfit12$CSWA <- Species$CSWA # r2 = 0.214
speciesfit12$PUMA <- Species$PUMA ##r2 = 0.3774
speciesfit12$COGR <- Species$COGR ##r2 = 0.2187
speciesfit12$AMRO <- Species$AMRO  ##r2 = 0.1589

speciesfit12


speciesfit12 <- envfit(species.nmds, speciesfit12, permutations = 999, choices = c(1,2))
speciesfit12

envfit12 <- as.data.frame(speciesfit12$vectors$arrows*sqrt(speciesfit12$vectors$r)) #scaling vectors
envfit12$factors<-rownames(speciesfit12)
envfit12


plot1 <-ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  theme_classic() +
  geom_point(shape = Year) + 
  stat_ellipse(data = data.scores,aes(x=NMDS1,y=NMDS2, linetype = vegYr), level = 0.90, size=1) +
  geom_segment(data=envfit12,aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length=unit(0.2,"cm"), type="closed"), colour="black") + #adding vectors
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  coord_equal() +
  theme(axis.title.x = element_text(margin=margin(t=10), size=18), 
        axis.title.y = element_text(margin=margin(r=10), size=18),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))
plot1

