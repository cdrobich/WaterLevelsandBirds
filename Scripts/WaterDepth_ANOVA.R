
##### Differences in Water Depth among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)
library(ggpubr)

Data <- read.csv("Data/2014_2015_univariate.csv")
str(Data)

colnames(Data)
unique(Data$Vegetation.type)


Univariate <- Data %>% #rename the factors
  mutate(Vegetation.type = fct_recode(Vegetation.type,
                                      "Emergent" = "Typha",
                                      "Invaded" = "Phragmites")) %>% 
  mutate(Year = fct_recode(Year,
                           "2014" = "Four",
                           "2015" = "Five"))

Univariate %>% count(Vegetation.type)


#### Two-way ANOVA comparing water depths ####

waterlevel <- lm(Water ~ Vegetation.type * Year, data = Univariate)

Anova(waterlevel, type = "3")

# Response: Water
# Sum Sq Df F value   Pr(>F)    
# (Intercept)           600.00  1  7.5500 0.009532 ** 
# Vegetation.type      1686.09  2 10.6083 0.000263 ***
# Year                  300.00  1  3.7750 0.060342 .  
# Vegetation.type:Year  297.66  2  1.8728 0.169211    
# Residuals            2701.99 34 


# Post-hoc test

hoc <- HSD.test(waterlevel, "Vegetation.type", group = TRUE, console = TRUE)

# Water groups
# Typha      20.14562      a
# Phragmites 17.25000      a
# Meadow      5.00000      b


plot(residuals(waterlevel)~fitted(waterlevel))


Univariate %>% group_by(Vegetation.type) %>% summarise(Water.avg = mean(Water),
                                                             Water.sd = sd(Water),
                                                             Water.min = min(Water),
                                                             Water.max = max(Water))

# Vegetation.type Year  Water.avg Water.sd Water.min Water.max
# 1 Meadow          Five      10        6.16         0        16
# 2 Meadow          Four       0        0            0         0
# 3 Phragmites      Five      28.2     13.9          5        42
# 4 Phragmites      Four       6.33     6.02         0        16
# 5 Typha           Five      31.1      8.15        21        42
# 6 Typha           Four       9.17    11.3          0        34


Univariate %>% group_by(Vegetation.type, Year) %>% summarise(Water.avg = mean(Water),
                                            Water.sd = sd(Water),
                                            Water.min = min(Water),
                                            Water.max = max(Water),
                                            range = (Water.max - Water.min))

#Vegetation.type    Year  Water.avg  Water.sd  Water.min  Water.max  range
# 1 Meadow         2015      10        6.16         0        16    16
#2 Meadow          2014       0        0            0         0     0
#3 Invaded         2015      28.2     13.9          5        42    37
#4 Invaded         2014       6.33     6.02         0        16    16
#5 Emergent        2015      31.1      8.15        21        42    21
#6 Emergent        2014       9.17    11.3          0        34    34

#Year  Water.avg Water.sd Water.min Water.max
# 2015     23.9     13.2          0        42
# 2014       5.57     8.49         0        34


## Water depth figures
Univariate$Year <- factor(Univariate$Year, levels = c("2014", "2015")) 

Water <- ggplot(Univariate, aes(x = Vegetation.type, y = Water))

Waters <- Water + geom_jitter(
  aes(shape = Year, color = Year), 
  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  size = 4) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_sdl", fun.args = list(mult = 1),
    geom = "pointrange", size = 0.6,
    position = position_dodge(0.8)
  ) +
  labs(x = " ",
       y = expression(paste("Water Depth (cm)")))


WaterDepth <- Waters + scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))

WaterDepth 

ggsave("Figures/Water Depth ANOVA.JPEG")

################ Lake Erie ######

erie <- read.csv("Data/lakeerie_2014_2015_1918_2013.csv")
erie
colnames(erie)
str(erie)




erie$Month <- factor(erie$Month, levels = c("May","June", "July ", "August"))

erie.figure <- ggplot(data = erie, aes(x = Month, y = Average, group = Year, 
                                       colour = Year, shape = Year)) +
  geom_line() +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = Average - CI, ymax = Average + CI),
                width = 0.3) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(fill = NA)) +
  ylim(173.9, 175.0) +
  xlab(" ") +
  ylab("Lake Erie Depth relative to IGLD 1985 (m) ") +
  theme(axis.text = element_text(size = 13),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 15)) +
  scale_color_manual(values = c("#bdbdbd","#fc8d62","#1f78b4"))

erie.figure

ggsave("Figures/Lake Erie.JPEG")


water.panel <- ggarrange(erie.figure, WaterDepth, 
                         legend = "bottom",
                         widths = c(1,1),
                         heights = c(1,1),
                         align = "h")
                         
water.panel

ggsave("Figures/Sites_Erie.JPEG")
