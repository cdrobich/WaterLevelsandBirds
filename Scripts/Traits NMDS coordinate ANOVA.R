
##### Differences in richness and abundance among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)
library(gridExtra) # arranging panels

#### Traits NMDS coordinates ####

Traits <- read.csv("Data/traits_NMDS coordiantes.csv")
glimpse(Traits)

colnames(Traits)

Traits$Year <- as.factor(Traits$Year)

unique(Traits$VegType)


#### Histograms ####

par(mfrow = c(2,2))

hist(Traits$Axis1, 
     xlab = "Axis1 scores", main = " ",
     border = "black",
     col = "white")

hist(Traits$Axis1, 
     xlab = "Traits$Axis2", main = " ",
     border = "black",
     col = "white")

hist(Traits$Axis1, 
     xlab = "Traits$Axis3", main = " ",
     border = "black",
     col = "white")


#### Traits Axis 1 coordinate ANOVA ####
colnames(Traits)

Axis1ANOVA <- lm(Axis1 ~ VegYr, data = Traits)

Anova(Axis1ANOVA, type = "3")

# Response: Axis1
#             Sum Sq Df F value    Pr(>F)    
# (Intercept) 0.3717  1  1.3271 0.2573502    
# VegYr       7.8258  5  5.5891 0.0007308 ***
# Residuals   9.5213 34                      


hoc <- HSD.test(Axis1ANOVA, "VegYr", group = TRUE, console = TRUE)    


#Axis1 groups

# Meadow_2015    0.8980333      a
# Invaded_2014   0.1198500     ab
# Emergent_2015  0.0934875     ab
# Emergent_2014 -0.2155375      b
# Invaded_2015  -0.2788667      b
# Meadow_2014   -0.5762833      b


plot(residuals(Axis1ANOVA)~fitted(Axis1ANOVA))


## Figure
colnames(Traits)

Axis1Fig <- ggplot(Traits, aes(x = VegType, y = Axis1)) +
  geom_jitter(aes(shape = Year, color = Year), 
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
       y = expression(paste("NMDS Axis 1 Coordinates"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")

Axis1Fig





#### Traits Axis 1 coordinate ANOVA ####
colnames(Traits)

Axis2ANOVA <- lm(Axis2 ~ VegYr, data = Traits)

Anova(Axis2ANOVA, type = "3")

# Response: Axis2
# Sum Sq Df F value  Pr(>F)  
# (Intercept)  0.5030  1  1.5422 0.22280  
# VegYr        3.8921  5  2.3868 0.05852 .
# Residuals   11.0886 34

plot(residuals(Axis2ANOVA)~fitted(Axis2ANOVA))


## Figure
colnames(Traits)

Axis2Fig <- ggplot(Traits, aes(x = VegType, y = Axis2)) +
  geom_jitter(aes(shape = Year, color = Year), 
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
       y = expression(paste("NMDS Axis 2 Coordinates"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")

Axis2Fig


##### NMDS Axis 3 ####

colnames(Traits)

Axis3ANOVA <- lm(Axis3 ~ VegYr, data = Traits)

Anova(Axis3ANOVA, type = "3")

#Response: Axis3
#Sum Sq Df F value  Pr(>F)  
#(Intercept) 0.8665  1  4.7592 0.03615 *
# VegYr       1.4815  5  1.6274 0.17934  
# Residuals   6.1905 34 

           

plot(residuals(Axis3ANOVA)~fitted(Axis3ANOVA))


## Figure
colnames(Traits)

Axis3Fig <- ggplot(Traits, aes(x = VegType, y = Axis3)) +
  geom_jitter(aes(shape = Year, color = Year), 
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
       y = expression(paste("NMDS Axis 3 Coordinates"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")

Axis3Fig


### panels

grid.arrange(Axis1Fig, Axis2Fig, Axis3Fig, ncol = 3)

coordinate <- arrangeGrob(Axis1Fig, Axis2Fig, Axis3Fig, ncol = 3)

ggsave("Figures/Traits NMDS coordinate ANOVAS.jpeg", coordinate)
