library(vegan)

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

ggplot(species, aes(x = VegType, y = Shannon)) +
  geom_boxplot(aes(fill = Year), 
              position = position_dodge(0.9)) +
   theme_classic() +
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(x = " ",
       y = expression(paste("Shannon Diversity Index (H)")))
  
  

ggplot(species, aes(x = VegType, y = Pielous)) +
  geom_boxplot(aes(fill = Year), 
               position = position_dodge(0.9)) +
  theme_classic() +
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(x = " ",
       y = expression(paste("Pielou Evenness (J)"))) +
  ylim(0.5, 1.2)

