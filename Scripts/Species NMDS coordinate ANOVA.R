
##### Differences in richness and abundance among sites between years #####

library(vegan)
library(agricolae) #skewness, kurtosis, Tukeys
library(tidyverse)
library(car) ## ANOVA function
library(Hmisc)
library(gridExtra) # arranging panels

#### Species NMDS coordinates ####

Species <- read.csv("Data/NMDS coordinates_species.csv")
glimpse(Species)

colnames(Species)

Species$Year <- as.factor(Species$Year)

unique(Species$VegType)


#### Histograms ####

par(mfrow = c(2,2))

hist(Species$Axis1, 
     xlab = "Axis1 scores", main = " ",
     border = "black",
     col = "white")

hist(Species$Axis1, 
     xlab = "Species$Axis2", main = " ",
     border = "black",
     col = "white")

hist(Species$Axis1, 
     xlab = "Species$Axis3", main = " ",
     border = "black",
     col = "white")


#### Species Axis 1 coordinate ANOVA ####
colnames(Species)

Axis1ANOVA <- lm(Axis1 ~ VegYr, data = Species)

Anova(Axis1ANOVA, type = "3")

#Response: Axis1
#            Sum Sq Df F value   Pr(>F)   
#(Intercept) 1.1716  1  4.0643 0.051756 . 
#VegYr       6.0795  5  4.2180 0.004329 **
# Residuals   9.8011 34 

hoc <- HSD.test(Axis1ANOVA, "VegYr", group = TRUE, console = TRUE)    

#Axis1 groups
#Emergent_2014  0.3826875      a
#Emergent_2015  0.2712500      a
#Invaded_2015   0.1734000      a
#Invaded_2014  -0.1096333     ab
#Meadow_2014   -0.1214167     ab
#Meadow_2015   -0.8142833      b


plot(residuals(Axis1ANOVA)~fitted(Axis1ANOVA))


## Figure
colnames(Species)

Axis1Fig <- ggplot(Species, aes(x = VegType, y = Axis1)) +
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
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(-2, 2)

Axis1Fig





#### Species Axis 1 coordinate ANOVA ####
colnames(Species)

Axis2ANOVA <- lm(Axis2 ~ VegYr, data = Species)

Anova(Axis2ANOVA, type = "3")

# Response: Axis2
#             Sum Sq Df F value  Pr(>F)   
# (Intercept) 0.0157  1  0.0653 0.79980   
# VegYr       4.5961  5  3.8207 0.00746 **
# Residuals   8.1800 34   

hoc <- HSD.test(Axis2ANOVA, "VegYr", group = TRUE, console = TRUE)

# Axis2 groups
# Meadow_2015    0.39788333      a
# Emergent_2015  0.37268750      a
# Emergent_2014 -0.04432500     ab
# Invaded_2015  -0.05896667     ab
# Invaded_2014  -0.14893333     ab
# Meadow_2014   -0.62788333      b

plot(residuals(Axis2ANOVA)~fitted(Axis2ANOVA))


## Figure
colnames(Species)

Axis2Fig <- ggplot(Species, aes(x = VegType, y = Axis2)) +
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
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") + 
  ylim(-2, 2)

Axis2Fig


##### NMDS Axis 3 ####

colnames(Species)

Axis3ANOVA <- lm(Axis3 ~ VegYr, data = Species)


Anova(Axis3ANOVA, type = "3")

hoc <- HSD.test(Axis2ANOVA, "VegYr", group = TRUE, console = TRUE)

#Response: Axis3
#Sum Sq Df F value Pr(>F)
#(Intercept)  0.0371  1  0.1175 0.7338
#VegYr        0.6212  5  0.3939 0.8495
#Residuals   10.7218 34                    


plot(residuals(Axis3ANOVA)~fitted(Axis3ANOVA))


## Figure
colnames(Species)

Axis3Fig <- ggplot(Species, aes(x = VegType, y = Axis3)) +
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
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(-2, 2)

Axis3Fig


### panels

grid.arrange(Axis1Fig, Axis2Fig, Axis3Fig, legend, widths = c(2.3, 2.3, 2.3, 0.8))

coordinate <- arrangeGrob(Axis1Fig, Axis2Fig, Axis3Fig, legend, widths = c(2.3, 2.3, 2.3, 0.8))

ggsave("Figures/Species NMDS coordinate ANOVAS.jpeg", coordinate)

### get the legend 

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


legend <- get_legend(Axis3Fig)
