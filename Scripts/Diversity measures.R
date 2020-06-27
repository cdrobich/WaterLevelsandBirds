library(vegan)
library(tidyverse)
library(gridExtra)

#### beta diversity  ####

env <- read.csv("Data/Env_variables.csv")
species <- read.csv("Data/Species matrix_raw.csv")

spp <- species[ , 5:27]

beta <- betadiver(spp, "-1")
mod <- with(env, betadisper(beta, VegYr))

anova(mod)

boxplot(mod,
        xlab = " ",
        ylab = "Distance to Centroid (Beta Diversity)")
       


#### alpha diversity ####
species$Year <- as.factor(species$Year)
str(species)

species$Shannon <- diversity(spp) # Shannon index

species$Pielous <- H/log(specnumber(spp)) # Pielou's evenness

hist(species$Shannon)
hist(species$Pielous)

Div <- ggplot(species, aes(x = VegType, y = Shannon)) +
  geom_boxplot(aes(fill = Year), 
              position = position_dodge(0.9)) +
   theme_classic() +
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(x = " ",
       y = expression(paste("Shannon Diversity Index (H)")))  +
  theme(legend.position = "blank") +
  ylim(0, 2)
  
  

Even <- ggplot(species, aes(x = VegType, y = Pielous)) +
  geom_boxplot(aes(fill = Year), 
               position = position_dodge(0.9)) +
  theme_classic() +
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(x = " ",
       y = expression(paste("Pielou Evenness (J)"))) +
  ylim(0, 2) +
  theme(legend.position = "blank")
Even
Div



grid.arrange(Div, Even, legend, ncol = 3, widths = c(2.3, 2.3, 0.8))

alphadiversity <- arrangeGrob(Div, Even, legend, ncol = 3, widths = c(2.3, 2.3, 0.8))

ggsave("Figures/Alphda diversity panel.jpeg", alphadiversity)


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(Even)
