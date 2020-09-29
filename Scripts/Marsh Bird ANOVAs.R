
##### Differences in richness and abundance among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)
library(gridExtra)
library(ggpubr)

Data <- read.csv("Data/2014_2015_univariate.csv")
str(Data)

dim(Data)
colnames(Data)
unique(Data$Vegetation.type)


Univariate <- Data %>% #rename the factors
  mutate(Vegetation.type = fct_recode(Vegetation.type,
                                      "Emergent" = "Typha",
                                      "Invaded" = "Phragmites")) 

Univariate %>% count(Vegetation.type, Year)


Univariate <- Univariate %>% #rename the factors
  mutate(Year = fct_recode(Year,
                           "2014" = "Four",
                           "2015" = "Five"))

Univariate$Year <- factor(Univariate$Year, levels = c("2014", "2015")) 

Univariate %>% count(Year)


colnames(Univariate)
str(Univariate)

Transform <- Univariate %>% mutate(logTS = log(TS),
                                   logTab= log(Tab),
                                   logMAb = log(Mab),
                                   logMS = log(MS))

Transform <- Transform %>% unite("VegYr", Vegetation.type:Year, remove = FALSE)


#### Marsh Abundance ANOVA w/log transformed & raw data ####

MabANOVA <- lm(logMAb ~ Vegetation.type * Year, data = Transform)

Anova(MabANOVA, type = "3")

# Response: logMAb
#                      Sum Sq  Df  F value    Pr(>F)    
# (Intercept)          17.3793  1  62.4974  3.307e-09 ***
# Vegetation.type       1.3266  2  2.3853    0.1073    
# Year                  0.6633  1  2.3854    0.1317    
# Vegetation.type:Year  0.3659  2  0.6579    0.5244    
# Residuals             9.4548 34                      


plot(residuals(MabANOVA)~fitted(MabANOVA))

## Total Marsh Abundance, raw data

MabANOVA2 <- lm(Mab ~ Vegetation.type * Year, data = Transform)


Anova(MabANOVA2, type = "3")

#Response: Mab
#                      Sum Sq Df F value    Pr(>F)    
# (Intercept)          280.17  1 25.4584 1.499e-05 ***
# Vegetation.type       52.78  2  2.3982    0.1061    
# Year                  24.08  1  2.1884    0.1483    
# Vegetation.type:Year  15.19  2  0.6902    0.5084    
# Residuals            374.17 34                        


Transform %>% group_by(Vegetation.type, Year) %>% summarise(MarshAb.avg = mean(Mab),
                                                            MarshAb.sd = sd(Mab),
                                                            MarshAb.min = min(Mab),
                                                            MarshAb.max = max(Mab))

#Vegetation.type Year  MarshAb.avg MarshAb.sd MarshAb.min MarshAb.max
#1 Meadow          Five         6.83       4.96           3          14
#2 Meadow          Four         4          2              1           6
#3 Invaded         Five         7.67       3.33           3          13
#4 Invaded         Four         8          3.90           2          13
#5 Emergent        Five        10.5        2.98           6          15
#6 Emergent        Four         9.5        2.33           7          14



#### Marsh Abundance Richness figures ####

MarshAb <- ggplot(Transform, aes(x = Vegetation.type, y = Mab))

MarshAbundance <- MarshAb + geom_jitter(
  aes(shape = Year, color = Year), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 2) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", size = 0.6,
    position = position_dodge(0.8)) +
  labs(x = " ",
       y = expression(paste("Marsh Bird Abundance"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  ylim(0, 30) +
  theme(legend.position = "blank")


MarshAbundance




#### Marsh Species Richness ANOVA w/ raw data ####

MSANOVA <- lm(logMS ~ Vegetation.type * Year, data = Transform)

Anova(MSANOVA, type = "3")

#Response: logMS
#                      Sum Sq Df  F value  Pr(>F)    
#(Intercept)          3.9237  1  34.4911  1.259e-06 ***
# Vegetation.type     0.3403  2  1.4958   0.23844    
#Year                 0.3603  1  3.1675   0.08406 .  
#Vegetation.type:Year 0.2463  2  1.0826   0.35011    
#Residuals            3.8678 34


## Total Marsh Abundance, raw data

MSANOVA2 <- lm(MS ~ Vegetation.type * Year, data = Transform)


Anova(MSANOVA2, type = "3")

#Response: MS
#                     Sum Sq Df  F value    Pr(>F)    
#(Intercept)          37.500  1  50.1639   3.512e-08 ***
#Vegetation.type       1.342  2  0.8974    0.4171    
#Year                  2.083  1  2.7869    0.1042    
#Vegetation.type:Year  1.767  2  1.1816    0.3191    
#Residuals            25.417 34   

plot(residuals(MSANOVA2)~fitted(MSANOVA2))



Transform %>% group_by(Vegetation.type, Year) %>% summarise(MarshS.avg = mean(MS),
                                                            MarshS.sd = sd(MS),
                                                            MarshS.min = min(MS),
                                                            MarshS.max = max(MS))

#Vegetation.type Year  MarshS.avg MarshS.sd MarshS.min MarshS.max

#1 Meadow          Five        2.5      1.22           1          4
#2 Meadow          Four        1.67     0.516          1          2
#3 Invaded         Five        2.83     0.408          2          3
#4 Invaded         Four        3        0.894          2          4
#5 Emergent        Five        3.12     0.641          2          4
#6 Emergent        Four        3.12     1.13           2          5



#### Marsh Species Richness figures ####

MarshS <- ggplot(Transform, aes(x = Vegetation.type, y = MS))

MarshRichness <- MarshS + geom_jitter(
  aes(shape = Year, color = Year), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 2) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", size = 0.6,
    position = position_dodge(0.8)) +
  labs(x = " ",
       y = expression(paste("Marsh Species Richness"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  ylim(0, 10) 

MarshRichness


#### Panel ####

grid.arrange(TotalAbundance, TotalRichness,
             ncol = 2)

panel <- arrangeGrob(TotalAbundance, TotalRichness, MarshAbundance, MarshRichness)


ggsave("Figures/BirdUnivariate_panels.jpeg", panel)