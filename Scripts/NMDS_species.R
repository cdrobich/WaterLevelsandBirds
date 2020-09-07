library(vegan)
library(tidyverse)
library(ggpubr)

species <- read.csv("Data/Species matrix_column relativized.csv") # put them in the same order
species <- species[order(species$Site),]

dim(species)

(spp <- species[ , 7:29])
sp.env <- species[, 1:6]

sp.env$Year <- as.factor(sp.env$Year)
sp.env$VegYr <- as.factor(sp.env$VegYr)

str(sp.env)

sp.env.uninv <- sp.env %>% #rename the factors
  mutate(VegType = fct_recode(VegType,
                              "Uninvaded" = "Emergent",
                              "Uninvaded" = "Meadow")) 

sp.env.uninv <- sp.env.uninv %>% unite("Veg.Year", Year:VegType, remove = FALSE)

#### NMDS Ordination ####


k_vec <- 1:10 #dim 1 - 10

stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(spp, trace = FALSE)
set.seed(10)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i,
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3 axes looks best

### NMDS ordination ####

set.seed(126) 

spp.nms <- metaMDS(spp, distance = "bray",
                   autotransform = FALSE,
                   k = 3, trymax = 1000) 

spp.nms

#global Multidimensional Scaling using monoMDS

#Data:     spp 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1555859 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘spp’

spp.nms$iters # 192
spp.nms$ities
spp.nms$stress

# quick look at the NMDS and the Shepard plot (for fit)

layout(matrix(1:2, ncol = 2))
plot(spp.nms, main = "Bird Spp NMDS plot"); stressplot(spp.nms, main = "Shepard plot")
layout(1)

# non-metric fit, r2 = 0.976
# linear fit, r2 = 0.851

# Goodness of fit

(g <- goodness(spp.nms)) # smaller the number the better the fit
sum(g^2)
spp.nms$stress^2 # 0.0242


# Extract the scores

bird.scores <- as.data.frame(scores(spp.nms, display = "sites"))

bird.scores$Sites <- sp.env.uninv$Site
bird.scores$Year <- as.factor(sp.env.uninv$Year)
bird.scores$Vegetation <- sp.env.uninv$VegType
bird.scores$VegYr <- sp.env.uninv$Veg.Year

colnames(bird.scores)
str(bird.scores)

write.csv(bird.scores, "Data/NMDS_scores_birdspecies.csv")

## Vectors

set.seed(126)

#### put meadow and cattail together

unique(sp.env.uninv$VegType)
spp

### vectors correlated with axis 1 & 2
vector.12 <- envfit(spp.nms, spp, 
             choices = 1:2,
             permutation = 999)

vector.12
vector12.df <- data.frame((vector.12$vectors)$arrows, (vector.12$vectors)$r, (vector.12$vectors)$pvals)

write.csv(vector12.df, "Data/NMDS_species_vectors_12.csv")

vector.12$vectors$r[vector.12$vectors$r > 0.25] # r2 over 0.25

vector12.df

corr.sp.12 <- spp %>% select(MAWR,VIRA, RWBL, YWAR, SOSP, TRES, EAKI, PUMA)

corrtaxa.12 <- envfit(spp.nms$points, corr.sp.12,
                   permutations = 999, choices = c(1,2))


corrtaxa.12

# make a new data frame for the figure
corr.spp.12 <- as.data.frame(corrtaxa.12$vectors$arrows*sqrt(corrtaxa.12$vectors$r)) #scaling vectors so they correspond with r2
corr.spp.12$species <- rownames(corr.spp.12)



## Axis 1 and 3 vectors

vector.13 <- envfit(spp.nms, spp, 
                    choices = 1:3,
                    permutation = 999)

vector.13
vector13.df <- data.frame((vector.13$vectors)$arrows, (vector.13$vectors)$r, (vector.13$vectors)$pvals)

write.csv(vector13.df, "Data/NMDS_species_vectors_13.csv")

vector.13$vectors$r[vector.13$vectors$r > 0.25] # r2 over 0.25

corr.sp.13 <- species %>% select(MAWR, SORA, SWSP, VIRA, RWBL, YWAR, SOSP, TRES, BARS, EAKI, PUMA)

corrtaxa.13 <- envfit(spp.nms$points, corr.sp.13,
                      permutations = 999, choices = c(1,3))


corrtaxa.13

# make a new data frame for the figure
corr.spp.13 <- as.data.frame(corrtaxa.13$vectors$arrows*sqrt(corrtaxa.13$vectors$r)) #scaling vectors so they correspond with r2
corr.spp.13$species <- rownames(corr.spp.13)

#### Figures #######

bird.scores # scores
corr.spp.12 #vectors 1, 2
corr.spp.13 # vectors axis 1, 3
colnames(bird.scores)

## Axis 1, 2 

birdspp.12 <- ggplot(data = bird.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = bird.scores, aes(x = NMDS1, y = NMDS2, colour = Vegetation, shape = VegYr), size = 5) + # sites as points
  stat_ellipse(data = bird.scores, aes(x = NMDS1,y = NMDS2,linetype = VegYr, colour = Vegetation), size = 1) + # a 95% CI ellipses
  geom_segment(data = corr.spp.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_classic(base_size = 16) + # no background
  scale_color_manual(values = c("#01665e", "#8c510a")) + # adding colours
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") + 
  ylim(-2, 2) +
  geom_label(data = corr.spp.12, aes(x = MDS1,y = MDS2,label = species),size=4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  scale_linetype_manual(values = c(1, 1, 2, 2))

birdspp.12


## NMDS Axis 1, 3
# same as above
birdspp.13 <- ggplot(data = bird.scores,
                     aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = bird.scores, aes(x = NMDS1, y = NMDS3, colour = Vegetation, shape = VegYr), size = 5) + # sites as points
  stat_ellipse(data = bird.scores, aes(x = NMDS1,y = NMDS3,linetype = VegYr, colour = Vegetation), size = 1) + # a 95% CI ellipses
  geom_segment(data = corr.spp.13, aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_classic(base_size = 16) + # no background
  scale_color_manual(values = c("#01665e", "#8c510a")) + # adding colours
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") + 
  ylim(-2.5, 2.5) +
  geom_label(data = corr.spp.13, aes(x = MDS1,y = MDS3,label = species),size=5) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  scale_linetype_manual(values = c(1, 1, 2, 2))
  

birdspp.13


## this is from ggpubr
# arranging the two figures
NMDS.spp <- ggarrange(birdspp.12, birdspp.13, # put the plot items in
                      nrow = 2, # I want them on top of each other
                      common.legend = TRUE, # they have the same legend
                      legend = "none")

NMDS.spp

ggsave("Figures/NMDS_spps.jpeg", NMDS.spp,
       width = 10, height = 10, dpi = 150, units = "in") # save that figure to my folder

ggsave("Figures/NMDS_spps.jpeg", NMDS.spp2,
       width = 10, height = 10, dpi = 150, units = "in")



## beta disper

spp.b <- vegdist(spp, method = "bray")
groups <- sp.env.uninv$Veg.Year

(betadisper(spp.b, groups))
anova(betadisper(spp.b, groups))

plot(betadisper(spp.b, groups))
boxplot(betadisper(spp.b, groups))








