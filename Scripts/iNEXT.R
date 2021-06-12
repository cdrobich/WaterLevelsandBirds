
# library -----------------------------------------------------------------

library(iNEXT)
library(tidyverse)
library(patchwork)


# data --------------------------------------------------------------------
richness <- read.csv("Data/species_richness.csv")

remnant <- richness %>% 
  select(remnant2014, remnant2015)

invaded <- richness %>% 
  select(invaded2014, invaded2015)
invaded <- invaded[1:6,]

# remnant vegetation

x <- iNEXT(remnant, q = 0, datatype = "abundance")
y <- iNEXT(remnant, q = 1, datatype = "abundance")
z <- iNEXT(remnant, q = 2, datatype = "abundance")


#$AsyEst: asymptotic diversity estimates along with related statistics.
#Site         Diversity Observed Estimator  s.e.    LCL    UCL
#1 remnant2014  Species richness   14.000    14.000 0.419 14.000 15.024
#2 remnant2014 Shannon diversity   13.655    15.127 0.690 13.774 16.479
#3 remnant2014 Simpson diversity   13.280    16.374 1.316 13.795 18.953
#4 remnant2015  Species richness   14.000    14.000 0.264 14.000 14.592
#5 remnant2015 Shannon diversity   13.483    14.549 0.603 13.483 15.732
#6 remnant2015 Simpson diversity   12.993    15.071 0.897 13.313 16.829



x.plot <- ggiNEXT(x, type = 1, se = TRUE)+
  labs(y = 'Species Richness') +
  ggtitle("Remnant Vegetation")

y.plot <- ggiNEXT(y, type = 1, se = TRUE)+
  labs(y = 'Shannon Diversity') 

z.plot <- ggiNEXT(z, type = 1, se = TRUE) +
  labs(y = "Simpson's Diversity") 

rem.plot <- x.plot + y.plot + z.plot

# invaded vegetation

a <- iNEXT(invaded, q = 0, datatype = "abundance")
b <- iNEXT(invaded, q = 1, datatype = "abundance")
c <- iNEXT(invaded, q = 2, datatype = "abundance")


#$AsyEst: asymptotic diversity estimates along with related statistics.
#Site         Diversity Observed Estimator  s.e.   LCL   UCL
#1 invaded2014  Species richness    6.000     6.000 0.245 6.000 6.550
#2 invaded2014 Shannon diversity    5.725     6.165 0.387 5.725 6.924
#3 invaded2014 Simpson diversity    5.493     6.330 0.671 5.493 7.646
#4 invaded2015  Species richness    6.000     6.000 0.220 6.000 6.486
#5 invaded2015 Shannon diversity    5.908     6.424 0.297 5.908 7.007
#6 invaded2015 Simpson diversity    5.824     6.940 0.633 5.824 8.180


a.plot <- ggiNEXT(a, type = 1, se = TRUE) +
  labs(y = 'Species Richness') +
  ggtitle("Invaded vegetation") 

b.plot <- ggiNEXT(b, type = 1, se = TRUE) +
  labs(y = 'Shannon Diversity') 

c.plot <- ggiNEXT(c, type = 1, se = TRUE) +
  labs(y = "Simpson's Diversity") 

inv.plot <- a.plot + b.plot + c.plot


