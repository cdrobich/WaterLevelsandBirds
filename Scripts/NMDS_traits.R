library(vegan)
library(tidyverse)

traits <- read.csv("Data/traits_matrix_Routput.csv")
traits <- traits[ ,2:21]
env <- read.csv("Data/Env_variables.csv")

env <- env %>% unite("VegYr", Year:VegType, remove = FALSE)


env$Year <- as.factor(sp.env$Year)
env$VegYr <- as.factor(sp.env$VegYr)

str(sp.env)

#### NMDS Ordination ####

set.seed(10) 
trt.nms <- metaMDS(traits) 

(trt.nms2 <- metaMDS(traits, previous.best = trt.nms, trymax = TRUE)) 
Call:
  metaMDS(comm = traits, trymax = TRUE, previous.best = trt.nms) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(traits) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1648839 
#Stress type 1, weak ties
#Two convergent solutions found after 40 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(traits)’ 


#####################################################

#Scaling: centring, PC rotation, halfchange scaling 1
#Species: expanded scores based on ‘wisconsin(spp)’ 1

# quick look at the NMDS and the Shepard plot (for fit)
layout(matrix(1:2, ncol = 2))
plot(trt.nms2, main = "Bird Trait NMDS plot"); stressplot(trt.nms2, main = "Shepard plot")
layout(1)

# non-metric fit, r2 = 0.973
# linear fit, r2 = 0.9

#Plot stress of each dimension
k_vec <- 1:10 #dim 1 - 10

stress <- numeric(length(k_vec)) # stress of each model put here

dune_dij <- metaMDSdist(traits, trace = FALSE)

set.seed(10)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i,
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # lookin good after 2 dimensions


# Goodness of fit

(g <- goodness(trt.nms2)) # smaller the number the better the fit
summary(g)

set.seed(10)

ev <- envfit(trt.nms2 ~ VegYr, data = env, # the . means give me everything in the data
             choices = 1:2,
             scaling = 3, # could change to sites or species
             permutation = 1000)

ev

#**FACTORS:
  
#  Centroids:
#                   NMDS1   NMDS2
#VegYr2014_Emergent  0.1932  0.1035
#VegYr2014_Invaded  -0.1149  0.1039
#VegYr2014_Meadow    0.0590 -0.2097
#VegYr2015_Emergent  0.0860  0.0829
#VegYr2015_Invaded   0.1097 -0.1087
#VegYr2015_Meadow   -0.4261 -0.0341

#Goodness of fit:
#         r2  Pr(>r)  
#VegYr 0.2732 0.01299 *
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 1000


## GAM smoothing for water depth

surf <- ordisurf(trt.nms2 ~ Depth,
                 data = env,
                 knots = 10, # fairly flexible surface but only as smooth as it needs to be
                 isotropic = TRUE, # same smoothing in both axes
                 main = NULL)

summary(surf)

#Link function: identity 

#Formula:
#  y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)

#Parametric coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   14.733      2.249   6.552 9.47e-08 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(x1,x2) 0.5545      9 0.085   0.263

#R-sq.(adj) =  0.0193   Deviance explained = 3.32%
#-REML = 161.04  Scale est. = 202.24    n = 40

#######################################

## Trying to plot


colvec <- c("black", "blue", "deeppink",
            "gray70", "aquamarine", "pink")

disp <- "sites"
scl <- 3 
shp <- c(1, 0, 2, 19, 15, 17)


cols <- with(env, colvec[VegYr])
shps <- with(env, shp[VegYr])

ordiplot(trt.nms2, type = "n",
         scaling = scl,
         display = disp)

points(trt.nms2)


lvl <- with(env, levels(VegYr))

legend("topright", legend = lvl,
       bty = "n", col = colvec,
       pch = shp) # create the legend by hand and map colours to levels

ordihull(trt.nms2, groups = env$VegYr, ## adds hulls around the most extreme points within each group
         col = colvec,
         scaling = scl,
         lwd = 2)

# species vectors

trt.fit <- envfit(trt.nms2, traits, permutations = 999)
head(trt.fit)

plot(trt.fit, p.max = 0.01, col = "black")

ordispider(spp.nms2, groups = sp.env$VegYr, # show labeled centroid & distance to that centroid within groups
           col = colvec,
           scaling = scl,
           label = TRUE)


ordiellipse(spp.nms2, groups = sp.env$VegYr, # "standard error of centroid" ellipse, can assess differences
            draw = "polygon", col = col_vec,
            scaling = scl, lwd = 2)


