
##### Differences in richness and abundance among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)
library(gridExtra)

Data <- read.csv("Data/2014_2015_univariate.csv")
str(Data)

colnames(Data)
unique(Data$Vegetation.type)


Univariate <- Data %>% #rename the factors
  mutate(Vegetation.type = fct_recode(Vegetation.type,
                                      "Emergent" = "Typha",
                                      "Invaded" = "Phragmites")) 

Univariate %>% count(Vegetation.type)


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

str(Transform)


par(mfrow = c(2,2))

hist(Univariate$logTab, 
     xlab = "Total abundance", main = " ",
     border = "black",
     col = "white")

hist(Univariate$logTS, 
     xlab = "Total Species Richness", main = " ",
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


Transform %>% group_by(Vegetation.type, Year) %>% summarise(TotalAb.avg = mean(Tab),
                                                            TotalAb.sd = sd(Tab),
                                                            TotalAb.med = median(Tab),
                                                            TotalAb.min = min(Tab),
                                                            TotalAb.max = max(Tab))

#Vegetation.type Year  TotalAb.avg TotalAb.sd TotalAb.min TotalAb.max

#1 Meadow          Five         23         6.23          13          31
#2 Meadow          Four         9.5       3.08           4          12
#3 Invaded         Five         14.8       2.23          11          17
#4 Invaded         Four         16.2       4.79          12          24
#5 Emergent        Five         22.4       5.24          15          32
#6 Emergent        Four         15.9       3.76          11          24
 

#### Total Abundance figures ####

TotalAb <- ggplot(Transform, aes(x = Vegetation.type, y = Tab))

TotalAbundance <- TotalAb + geom_jitter(
  aes(shape = Year, color = Year), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 2) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", size = 0.6,
    position = position_dodge(0.8)
  ) +
  labs(x = " ",
       y = expression(paste("Total Bird Abundance"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")


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
                                                            TotalS.max = max(TS))

#Vegetation.type Year  TotalS.avg TotalS.sd TotalS.min TotalS.max

#1 Meadow          Five        7.5      1.87           5         10
#2 Meadow          Four        4.5      0.837          3          5
#3 Invaded         Five        5.17     0.983          4          6
#4 Invaded         Four        5.83     1.94           3          9
#5 Emergent        Five        5.38     1.19           4          7
#6 Emergent        Four        4.88     1.36           4          8


#### Total Species Richness figures ####

TotalS<- ggplot(Transform, aes(x = Vegetation.type, y = TS))

TotalRichness <- TotalS + geom_jitter(
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
       y = expression(paste("Total Bird Species Richness"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")


TotalRichness



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



#### Total Species Richness figures ####

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

grid.arrange(TotalAbundance, TotalRichness, MarshAbundance, MarshRichness)

panel <- arrangeGrob(TotalAbundance, TotalRichness, MarshAbundance, MarshRichness)


ggsave("Figures/BirdUnivariate_panels.jpeg", panel)

### get the legend 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


legend <- get_legend(TotalRichness)
legend

