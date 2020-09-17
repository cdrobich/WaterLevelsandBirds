
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



#### Histograms ####

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

#### transform variables ####

Transform <- Univariate %>% mutate(logTS = log(TS),
                      logTab= log(Tab),
                      logMAb = log(Mab),
                      logMS = log(MS))

Transform <- Transform %>% unite("VegYr", Vegetation.type:Year, remove = FALSE)


str(Transform)


par(mfrow = c(2,2))

hist(Transform$logTab, 
     xlab = "Total abundance", main = " ",
     border = "black",
     col = "white")

hist(Transform$logTS, 
     xlab = "Total SpeciesRichness", main = " ",
     border = "black",
     col = "white")

hist(Univariate$logMab, 
     xlab = "Total marsh bird abundance", main = " ",
     border = "black",
     col = "white")

hist(Univariate$logMS, 
     xlab = "Total marsh bird abundance", main = " ",
     border = "black",
     col = "white")


#### Total Abundance ANOVA w/log transformed & raw data ####

TAbANOVA <- lm(logTab ~ Vegetation.type * Year, data = Transform)


Anova(TAbANOVA, type = "3")


#Response: logTab
#                      Sum Sq  Df  F value    Pr(>F)    
# (Intercept)          57.663  1  759.8577  < 2.2e-16 ***
#  Vegetation.type       0.688  2   4.5342  0.0179652 *  
#  Year                  2.488  1  32.7899  1.953e-06 ***
#  Vegetation.type:Year  1.438  2   9.4747  0.0005364 ***
#  Residuals             2.580 34                       


plot(residuals(TAbANOVA)~fitted(TAbANOVA))

##

TAbANOVA2 <- lm(Tab ~ Vegetation.type * Year, data = Transform)


Anova(TAbANOVA2, type = "3")

# Response: Tab
#                       Sum Sq Df  F value    Pr(>F)    
#  (Intercept)          3174.0  1 160.6092 1.957e-14 ***
#  Vegetation.type       257.5  2   6.5147  0.004026 ** 
#  Year                  546.7  1  27.6664 7.899e-06 ***
#  Vegetation.type:Year  330.5  2   8.3608  0.001114 ** 
#  Residuals             671.9 34                       


Ab.sum <- Transform %>% group_by(Vegetation.type, Year) %>% summarise(TotalAb.avg = mean(Tab),
                                                            TotalAb.sd = sd(Tab),
                                                            N = length(VegYr),
                                                            TotalAb.se = (TotalAb.sd)/sqrt(N),
                                                            TotalAb.med = median(Tab),
                                                            TotalAb.min = min(Tab),
                                                            TotalAb.max = max(Tab))
Ab.sum <- as.data.frame(Ab.sum)
Ab.sum <- Ab.sum %>% unite("VegYr", Vegetation.type:Year, remove = FALSE)

#Vegetation.type Year  TotalAb.avg TotalAb.sd     N TotalAb.se TotalAb.med TotalAb.min TotalAb.max

#1 Meadow          2015         23         6.23     6      2.54         23.5          13          31
#2 Meadow          2014          9.5       3.08     6      1.26         10.5           4          12
#3 Invaded         2015         14.8       2.23     6      0.910        15            11          17
#4 Invaded         2014         16.2       4.79     6      1.96         14            12          24
#5 Emergent        2015         22.4       5.24     8      1.85         21.5          15          32
#6 Emergent        2014         15.9       3.76     8      1.33         15.5          11          24
 
Transform
#### Total Abundance figures ####


Transform <- Transform %>% 
  mutate(Year = fct_relevel(Year,
                            "2014",
                            "2015"))

Transform <- Transform %>% 
  mutate(Vegetation.type = fct_relevel(Vegetation.type,
                              "Invaded",
                              "Meadow",
                              "Emergent"))


TotalAbundance <- ggplot(Transform, aes(x = Vegetation.type, y = Tab)) + 
  geom_jitter(
    aes(shape = Year, color = Year), 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)
  ) +
  labs(x = " ",
       y = expression(paste("Total Bird Abundance"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 40)


TotalAbundance





#### Total Species Richness ANOVA w/log transformed & raw data ####

TSANOVA <- lm(logTS ~ Vegetation.type * Year, data = Transform)


Anova(TSANOVA, type = "3")

# Response: logTS
#                       Sum Sq   Df  F value    Pr(>F)    
#  (Intercept)          23.7063  1  378.2970 < 2.2e-16 ***
#  Vegetation.type       0.4950  2   3.9492  0.028694 *  
#  Year                  0.7519  1  11.9978  0.001459 ** 
#  Vegetation.type:Year  0.5437  2   4.3377  0.020996 *  
#  Residuals             2.1306 34                       


plot(residuals(TSANOVA)~fitted(TSANOVA))

## Total species richness, raw data

TSANOVA2 <- lm(TS ~ Vegetation.type * Year, data = Transform)


Anova(TSANOVA2, type = "3")

# Response: TS
#                       Sum Sq  Df  F value  Pr(>F)    
# (Intercept)          337.50   1 170.2101  8.59e-15 ***
#  Vegetation.type       20.74  2   5.2303  0.0104608 *  
#  Year                  27.00  1  13.6168  0.0007794 ***
#  Vegetation.type:Year  21.23  2   5.3543  0.0095171 ** 
#  Residuals             67.42 34                         


Transform %>% group_by(Vegetation.type, Year) %>% summarise(TotalS.avg = mean(TS),
                                                            TotalS.sd = sd(TS),
                                                            Total.med = median(TS),
                                                            TotalS.min = min(TS),
                                                            TotalS.max = max(TS),
                                                            N = length(VegYr),
                                                            TotalS.se = (TotalS.sd/sqrt(N)))

#Vegetation.type Year  TotalS.avg TotalS.sd Total.med TotalS.min TotalS.max     N   TotalS.se
#1 Meadow          2015        7.5      1.87        7.5          5         10     6     0.764
#2 Meadow          2014        4.5      0.837       5            3          5     6     0.342
#3 Invaded         2015        5.17     0.983       5.5          4          6     6     0.401
#4 Invaded         2014        5.83     1.94        6            3          9     6     0.792
#5 Emergent        2015        5.38     1.19        5            4          7     8     0.420
#6 Emergent        2014        4.88     1.36        4.5          4          8     8     0.479


#### Total Species Richness figures ####

TotalRichness <- ggplot(Transform, aes(x = Vegetation.type, y = TS)) +
  geom_jitter(aes(shape = Year, color = Year), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
  size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Total Bird Species Richness"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 15)


TotalRichness



ab.richn <- ggarrange(TotalAbundance, TotalRichness,
                         common.legend = TRUE, 
                         legend = "bottom",
                         widths = c(1,1),
                         heights = c(1,1),
                         align = "h",
                         labels = c("C","D"), 
                      hjust = c(-6, -6),
                      vjust = 2.5)
ab.richn

ggsave("Figures/Rich_Abun_panel.JPEG")


#### with remnant vs invaded

transform2 <- Transform %>% 
  mutate(Vegetation.type = fct_recode(Vegetation.type,
                                      "Remnant" = "Emergent",
                                      "Remnant" = "Meadow")) 


Abundance.remnant <- lm(Tab ~ Vegetation.type * Year, data = transform2)
Anova(Abundance.remnant, type = "3")

plot(Abundance.remnant$residuals)
plot(Abundance.remnant)


#                       Sum Sq Df F value    Pr(>F)    
#  (Intercept)          1568.17  1 69.4737 6.347e-10 ***
#  Vegetation.type        38.40  1  1.7013  0.200393    
#  Year                    5.33  1  0.2363  0.629852    
# Vegetation.type:Year  246.46  1 10.9187  0.002161 ** 
#  Residuals             812.60 36                      

Abundance.remnant2 <- lm(logTab ~ Vegetation.type * Year, data = transform2)
Anova(Abundance.remnant2, type = "3")
plot(Abundance.remnant2) # these look better

#                       Sum Sq Df  F value    Pr(>F)    
# (Intercept)          45.357  1 449.7785 < 2.2e-16 ***
#  Vegetation.type       0.250 1   2.4791  0.124117    
# Year                  0.012  1   0.1176  0.733698    
# Vegetation.type:Year  0.883  1   8.7581  0.005422 ** 
#  Residuals             3.630 36    

TotalAbundance2 <- ggplot(transform2, aes(x = Vegetation.type, y = Tab)) + 
  geom_jitter(
    aes(shape = Year, color = Year), 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)
  ) +
  labs(x = " ",
       y = expression(paste("Total Bird Abundance"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 40)


TotalAbundance2


Rich.remnant <- lm(TS ~ Vegetation.type * Year, data = transform2)
Anova(Rich.remnant, type = "3")
plot(Rich.remnant)

#                         Sum Sq Df F value    Pr(>F)    
#  (Intercept)          204.167  1 88.1496 3.258e-11 ***
#  Vegetation.type        5.260  1  2.2708   0.14056    
#  Year                   1.333  1  0.5757   0.45295    
# Vegetation.type:Year  10.519  1  4.5416   0.03998 *  
#  Residuals             83.381 36 


Rich.remnant2 <- lm(logTS ~ Vegetation.type * Year, data = transform2)
Anova(Rich.remnant2, type = "3")
plot(Rich.remnant2)

#Response: logTS
#                       Sum Sq Df  F value  Pr(>F)    
#(Intercept)          17.6150  1 252.2236 < 2e-16 ***
#  Vegetation.type     0.1462  1   2.0939 0.15654    
#Year                  0.0228  1   0.3267 0.57117    
#Vegetation.type:Year  0.2739  1   3.9220 0.05534 .  
#Residuals             2.5142 36                     


TotalRichness2 <- ggplot(transform2, aes(x = Vegetation.type, y = TS)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Total Bird Species Richness"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 15)


TotalRichness2



ab.richn2<- ggarrange(TotalAbundance2, TotalRichness2,
                      common.legend = TRUE, 
                      legend = "none",
                      widths = c(1,1),
                      heights = c(1,1),
                      align = "h",
                      labels = c("A","B"),
                      hjust = c(-6, -6),
                      vjust = 2.5)
ab.richn2

ggsave("Figures/Rich_Abun_panel_remnant_invaded.JPEG")




panel <- ggarrange(ab.richn2, ab.richn,
          nrow = 2)

ggsave("Figures/Rich_Abun_panel_two_three_groups.JPEG", panel)




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


########## Alpha Diversity Measures ##########

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
                              "Remnant" = "Emergent",
                              "Remnant" = "Meadow")) 

sp.env.uninv <- sp.env.uninv %>% unite("Veg.Year", Year:VegType, remove = FALSE)



### Diversity Measures

(H <- diversity(spp)) # Shannon Index
(J <- H/log(specnumber(spp))) # Pielou
(D <- diversity(spp, "simpson"))  # Simpsons 1 - D
(I <- diversity(spp, "invsimp")) # Simpsons 1/D
specnumber(spp) # species richness

diversity <- sp.env

diversity$H <- H
diversity$J <- J
diversity$D <- D
diversity$I <- I

diversity <- diversity %>% 
  mutate(VegType = fct_recode(VegType,
                              "Remnant" = "Emergent",
                              "Remnant" = "Meadow")) 

diversity <- diversity %>% unite("Veg.Year", Year:VegType, remove = FALSE)

diversity <- diversity %>% 
  mutate(VegType = fct_relevel(VegType,
                               "Invaded", "Remnant"))




diversity <- diversity %>% mutate(logSimp = log(I),
                                   logPielou= log(J))

par(mfrow = c(2,2))

hist(diversity$I, 
     xlab = "Simpsons", main = " ",
     border = "black",
     col = "white")

hist(diversity$logSimp, 
     xlab = "log Simpsons", main = " ",
     border = "black",
     col = "white")

hist(diversity$J, 
     xlab = "Pielou", main = " ",
     border = "black",
     col = "white")

hist(diversity$logPielou, 
     xlab = "log Pielou", main = " ",
     border = "black",
     col = "white")



shannon <- ggplot(diversity, aes(x = VegType, y = H)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Shannon-Weiner Diversity (H')"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 2)

shannon



simpson.anova <- lm(I ~ VegType * Year, data = diversity)
Anova(simpson.anova, type = "3")
plot(simpson.anova)

#Response: I
#              Sum Sq Df F value    Pr(>F)    
# (Intercept)  57.746  1 47.5141 4.537e-08 ***
# VegType       0.662  1  0.5443    0.4654    
# Year          0.846  1  0.6963    0.4095    
# VegType:Year  0.041  1  0.0334    0.8560    
# Residuals    43.752 36


simpson.anova2 <- lm(logSimp ~ VegType * Year, data = diversity)
Anova(simpson.anova2, type = "3")
plot(simpson.anova2) # looks worse


simpson <- ggplot(diversity, aes(x = VegType, y = I)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Inverse Simpsons (1/D)"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 6)

simpson # higher this value the higher the diversity



simpsons <- ggplot(diversity, aes(x = VegType, y = D)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Simpson's Diversity"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 1)

simpsons


pielou.anova <- lm(J ~ VegType * Year, data = diversity)
Anova(pielou.anova, type = "3")
plot(pielou.anova)

# Response: J
# Sum Sq Df F value    Pr(>F)    
# (Intercept)  3.2547  1 96.2758 1.029e-11 ***
# VegType      0.0031  1  0.0911    0.7645    
# Year         0.0630  1  1.8635    0.1807    
# VegType:Year 0.0157  1  0.4630    0.5006    
# Residuals    1.2170 36  

pielou <- ggplot(diversity, aes(x = VegType, y = J)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Pielou's Evenness (J)"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 1)

pielou



alpha <- ggarrange(TotalAbundance2, TotalRichness2,
                  simpson, pielou,
                      common.legend = TRUE, 
                      legend = "bottom",
                      widths = c(1,1),
                      heights = c(1,1),
                      align = "h",
                      labels = c("A","B", "C", "D"), 
                      hjust = c(-6, -6, -5, -7),
                      vjust = 2.5)
alpha

ggsave("Figures/alpha_diversity_abundance.TIFF", alpha)

##### with all three groups #########

diversity1 <- diversity %>% #rename the factors
  mutate(Vegetation.type = fct_recode(VegYr,
                                      "Invaded" = "Invaded_2014",
                                      "Invaded" = "Invaded_2015",
                                      "Meadow" = "Meadow_2014",
                                      "Meadow" = "Meadow_2015",
                                      "Emergent" = "Emergent_2014",
                                      "Emergent" = "Emergent_2015"))
                                     



diversity1 <- diversity1 %>% 
  mutate(VegYr = fct_relevel(Vegetation.type,
                                  "Invaded", "Meadow", "Emergent"))

## with all three groups ugh #

simpson2 <- ggplot(diversity1, aes(x = Vegetation.type, y = I)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Inverse Simpsons (1/D)"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 6)

simpson2 # higher this value the higher the diversity

pielou2 <- ggplot(diversity1, aes(x = Vegetation.type, y = J)) +
  geom_jitter(aes(shape = Year, color = Year), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Pielou's Evenness (J)"))) + 
  scale_color_manual(values = c("#fc8d62","#35978f")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 15)) +
  theme(legend.position = "blank") +
  ylim(0, 1)

pielou2




alpha.panel <- ggarrange(TotalAbundance, TotalRichness,
                         simpson2, pielou2, 
                      common.legend = TRUE, 
                      legend = "bottom",
                      widths = c(1,1),
                      heights = c(1,1),
                      align = "h",
                      labels = c("A", "B", "C","D"), 
                      hjust = c(-6, -6, -5, -6.5),
                      vjust = 2.5)
alpha.panel 

ggsave("Figures/alpha_panels.TIFF")
