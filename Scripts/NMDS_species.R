library(vegan)
library(tidyverse)

species <- read.csv("Data/Species matrix_raw.csv")

spp <- species[ , 5:27]
sp.env <- species[ , 1:4]

sp.env <- sp.env %>% unite("VegYr", Year:VegType, remove = FALSE)


sp.env$Year <- as.factor(sp.env$Year)
sp.env$VegYr <- as.factor(sp.env$VegYr)

str(sp.env)

#### NMDS Ordination ####

set.seed(10) 
spp.nms <- metaMDS(spp) 

(spp.nms2 <- metaMDS(spp, previous.best = spp.nms, trymax = TRUE)) 

#Call:
# metaMDS(comm = spp, trymax = TRUE, previous.best = spp.nms) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(spp) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.2002965   used default monoMDS so proportions on scale of 0-1

#Stress type 1, weak ties

#Two convergent solutions found after 40 tries

#Scaling: centring, PC rotation, halfchange scaling 1
#Species: expanded scores based on ‘wisconsin(spp)’ 1

# quick look at the NMDS and the Shepard plot (for fit)
layout(matrix(1:2, ncol = 2))
plot(spp.nms2, main = "Bird Spp NMDS plot"); stressplot(spp.nms2, main = "Shepard plot")
layout(1)

# non-metric fit, r2 = 0.96
# linear fit, r2 = 0.816

#Plot stress of each dimension
k_vec <- 1:10 #dim 1 - 10

stress <- numeric(length(k_vec)) # stress of each model put here

dune_dij <- metaMDSdist(spp, trace = FALSE)

set.seed(10)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i,
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # lookin good after 2 dimensions


# Goodness of fit

(g <- goodness(spp.nms2)) # smaller the number the better the fit
summary(g)

set.seed(10)

ev <- envfit(spp.nms2 ~ VegYr, data = sp.env, # the . means give me everything in the data
             choices = 1:2,
             scaling = "symmetric", # could change to sites or species
             permutation = 1000)

ev

#***FACTORS:
  
#  Centroids:
#                     NMDS1   NMDS2
#VegYr2014_Emergent  0.2347  0.1310
#VegYr2014_Invaded  -0.1831  0.0015
#VegYr2014_Meadow    0.1289 -0.4862
#VegYr2015_Emergent  0.1170  0.3312
#VegYr2015_Invaded   0.0406  0.0187
#VegYr2015_Meadow   -0.4554 -0.1501

# Goodness of fit:
#  r2  Pr(>r)  
# VegYr 0.2734 0.01299 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 1000




## GAM smoothing for water depth

surf <- ordisurf(spp.nms2 ~ Depth,
                 data = sp.env,
                 knots = 10, # fairly flexible surface but only as smooth as it needs to be
                 isotropic = TRUE, # same smoothing in both axes
                 main = NULL)

summary(surf)

#Family: gaussian 
#Link function: identity 

#Formula:
#  y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)

#Parametric coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   14.733      1.865   7.898 2.03e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F  p-value    
#s(x1,x2) 2.422      9 2.086 0.000324 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.325   Deviance explained = 36.7%
#-REML = 156.06  Scale est. = 139.2     n = 40



## Trying to plot


colvec <- c("black", "blue", "deeppink",
            "gray70", "aquamarine", "pink")

disp <- "sites"
scl <- 3 
shp <- c(1, 0, 2, 19, 15, 17)


cols <- with(sp.env, colvec[VegYr])
shps <- with(sp.env, shp[VegYr])

ordiplot(spp.nms2, type = "n",
     scaling = scl,
     display = disp)

points(spp.nms2, display = disp, scaling = scl,
       pch = shps, col = cols, cex = 2)


lvl <- with(sp.env, levels(VegYr))

legend("topright", legend = lvl,
       bty = "n", col = colvec,
       pch = shp) # create the legend by hand and map colours to levels

ordihull(spp.nms2, groups = sp.env$VegYr, ## adds hulls around the most extreme points within each group
         col = colvec,
         scaling = scl,
         lwd = 2)

# species vectors

spp.fit <- envfit(spp.nms2, spp, permutations = 999)
head(spp.fit)

plot(spp.fit, p.max = 0.01, col = "black")

ordispider(spp.nms2, groups = sp.env$VegYr, # show labeled centroid & distance to that centroid within groups
           col = colvec,
           scaling = scl,
           label = TRUE)


ordiellipse(spp.nms2, groups = sp.env$VegYr, # "standard error of centroid" ellipse, can assess differences
            draw = "polygon", col = col_vec,
            scaling = scl, lwd = 2)


