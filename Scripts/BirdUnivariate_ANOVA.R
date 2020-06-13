
##### Differences in richness and abundance among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)

Data <- read.csv("Data/2014_2015_univariate.csv")
str(Data)

colnames(Data)
unique(Data$Vegetation.type)


Univariate <- Data %>% #rename the factors
  mutate(Vegetation.type = fct_recode(Vegetation.type,
                                      "Emergent" = "Typha",
                                      "Invaded" = "Phragmites"))

Univariate %>% count(Vegetation.type)

colnames(Univariate)
str(Univariate)

Univariate$Mab <- as.integer(Univariate$Mab)

#### Histograms

par(mfrow = c(2,2))

hist(Univariate$Tab, 
     xlab = "Total abundance", main = " ",
     border = "black",
     col = "white")

hist(Univariate$TS, 
     xlab = "Total Species Richness", main = " ",
     border = "black",
     col = "white")


hist(Univariate$Mab, 
     xlab = "Total marsh bird abundance", main = " ",
     border = "black",
     col = "white")


hist(Univariate$MS, 
     xlab = "Total marsh bird richness", main = " ",
     border = "black",
     col = "white")

