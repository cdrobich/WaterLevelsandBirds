library(vegan)
library(tidyverse)

traits2 <- read.csv("Data/traits_matrix_rel_Routput.csv")
traits <- traits2[ ,2:21]
env <- read.csv("Data/Env_variables.csv")

env.t <- env %>% unite("VegYr", Year:VegType, remove = FALSE)


env.t$Year <- as.factor(env.t$Year)
env.t$VegYr <- as.factor(env.t$VegYr)

str(env)

# put uninvaded together into one group
trt.env.uninv <- env.t %>% #rename the factors
  mutate(VegType = fct_recode(VegType,
                              "Uninvaded" = "Emergent",
                              "Uninvaded" = "Meadow")) 

trt.env.uninv <- trt.env.uninv %>% unite("Veg.Year", Year:VegType, remove = FALSE)

#### NMDS Ordination ####

k_vec <- 1:10 #dim 1 - 10

stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(traits, trace = FALSE)
set.seed(10)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i,
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3 axes looks best

### NMDS ordination ####

set.seed(126) 

trt.nms <- metaMDS(traits, distance = "bray",
                   autotransform = FALSE,
                   k = 3, trymax = 1000) 

trt.nms

#Data:     traits 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.06959201 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘traits’


trt.nms$iters # 93


# quick look at the NMDS and the Shepard plot (for fit)

layout(matrix(1:2, ncol = 2))
plot(trt.nms, main = "Bird Spp NMDS plot"); stressplot(spp.nms, main = "Shepard plot")
layout(1)

# non-metric fit, r2 = 0.998
# linear fit, r2 = 0.993

# Goodness of fit

(g <- goodness(trt.nms)) # smaller the number the better the fit
sum(g^2)
trt.nms$stress^2 # 0.004


# Extract the scores

trt.scores <- as.data.frame(scores(trt.nms, display = "sites"))

trt.scores$Sites <- trt.env.uninv$Site
trt.scores$Year <- as.factor(trt.env.uninv$Year)
trt.scores$Vegetation <- trt.env.uninv$VegType
trt.scores$VegYr <- as.factor(trt.env.uninv$Veg.Year)

colnames(trt.scores)

write.csv(trt.scores, "Data/NMDS_scores_traits.csv")

## Vectors

set.seed(126)

#### put meadow and cattail together

unique(trt.env.uninv$VegType)

### vectors correlated with axis 1 & 2
vector.t.12 <- envfit(trt.nms, traits, 
                    choices = 1:2,
                    permutation = 999)

vector.t.12
vector12.t.df <- data.frame((vector.t.12$vectors)$arrows, (vector.t.12$vectors)$r, (vector.t.12$vectors)$pvals)

write.csv(vector12.t.df, "Data/NMDS_traits_vectors_12.csv")
vector.t.12$vectors$r[vector.t.12$vectors$r > 0.25] # r2 over 0.25


corr.t.12 <- traits %>% select(MARSH, OTHER, INSECT, FISH, 
                               GRD_FRG, STALK, AER_FORG, SHRUB, GROUND)

corrtraits.12 <- envfit(trt.nms$points, corr.t.12,
                      permutations = 999, choices = c(1,2))


# make a new data frame for the figure
corr.trt.12 <- as.data.frame(corrtraits.12$vectors$arrows*sqrt(corrtraits.12$vectors$r)) #scaling vectors so they correspond with r2
corr.trt.12$species <- rownames(corr.trt.12)



## Axis 1 and 3 vectors

vector.t.13 <- envfit(trt.nms, traits, 
                    choices = 1:3,
                    permutation = 999)

vector.t.13
vector13.t.df <- data.frame((vector.t.13$vectors)$arrows, (vector.t.13$vectors)$r, (vector.t.13$vectors)$pvals)
write.csv(vector13.t.df, "Data/NMDS_traits_vectors_13.csv")

vector.t.13$vectors$r[vector.t.13$vectors$r > 0.25] # r2 over 0.25

corr.trt.13 <- traits %>% select(MARSH, OTHER, INSECT, SEED, FISH, GRD_FRG, STALK, 
                                AER_FORG, SHRUB, GROUND, CAVITY, BUILD)

corrtraits.13 <- envfit(trt.nms$points, corr.trt.13,
                      permutations = 999, choices = c(1,3))


# make a new data frame for the figure
corr.trt.13 <- as.data.frame(corrtraits.13$vectors$arrows*sqrt(corrtraits.13$vectors$r)) #scaling vectors so they correspond with r2
corr.trt.13$species <- rownames(corr.trt.13)

#### Figures #######

trt.scores
corr.trt.12
corr.trt.13
str(corr.trt.12)
str(corr.trt.13)
## Axis 1, 2 

traits.12 <- ggplot(data = trt.scores,
                     aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = trt.scores, aes(x = NMDS1, y = NMDS2, colour = Vegetation, shape = VegYr), size = 5) + # sites as points
  stat_ellipse(data = trt.scores, aes(x = NMDS1,y = NMDS2,linetype = VegYr, colour = Vegetation), size = 1) + # a 95% CI ellipses
  geom_segment(data = corr.trt.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + 
  geom_label(data = corr.trt.12, aes(x = MDS1,y = MDS2,label = species),size=5) + # can add in geom_label or geom_text for labels
  theme_classic(base_size = 16) + # no background
  scale_color_manual(values = c("#01665e", "#8c510a")) + # adding colours
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") + 
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  scale_linetype_manual(values = c(1, 1, 2, 2))

traits.12

str(trt.scores)

## NMDS Axis 1, 3
traits.13 <- ggplot(data = trt.scores,
                     aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = trt.scores, aes(x = NMDS1, y = NMDS3, colour = Vegetation, shape = VegYr), size = 5) + # sites as points
  stat_ellipse(data = trt.scores, aes(x = NMDS1,y = NMDS3,linetype = VegYr, colour = Vegetation), size = 1) + # a 95% CI ellipses
  geom_segment(data = corr.trt.13, aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_classic(base_size = 16) + # no background
  scale_color_manual(values = c("#01665e", "#8c510a")) + # adding colours
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  ylim(-2, 2) +
  geom_label(data = corr.trt.13, aes(x = MDS1,y = MDS3,label = species),size=5) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  scale_linetype_manual(values = c(1, 1, 2, 2))

traits.13

## this is from ggpubr
# arranging the two figures
NMDS.trt <- ggarrange(traits.12, traits.13, # put the plot items in
                      nrow = 2, # I want them on top of each other
                      common.legend = TRUE, # they have the same legend
                      legend = "none")

NMDS.trt

ggsave("Figures/NMDS_traits.jpeg", NMDS.trt,
       width = 10, height = 10, dpi = 150, units = "in") # save that figure to my folder

NMDS.trt
NMDS.spp

panel <- ggarrange(NMDS.spp, NMDS.trt,
          legend = "none",
          labels = "AUTO",
          hjust = c(-6, -7),
          vjust = 2.5)
panel

ggsave("Figures/NMDS_panel.jpeg", panel,
       width = 11, height = 12, dpi = 150, units = "in") 