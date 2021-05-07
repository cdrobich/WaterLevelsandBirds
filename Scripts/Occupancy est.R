library(vegan)
library(tidyverse)


species <- read.csv("Data/Species matrix_raw.csv")

species$Year <- as.factor(species$Year)

species <- species %>% 
  mutate(Veg = fct_recode(VegType,
                          "Remnant" = "Emergent",
                          "Remnant" = "Meadow")) %>% 
  unite("VegYr", Year:VegType, remove = FALSE)

spp <- species %>% select(AMBI:AMRO)


# ChaoSpecies -------------------------------------------------------------

install.packages("SpadeR")
library(SpadeR)

Diversity(spp, "abundance")
