
library(vegan)
library(tidyverse)
library(ggpubr)
library(adespatial)
library(betapart)
library(picante) # randomizing matrix

### Import Bird matrix ####

species <- read.csv("Data/Species matrix_raw.csv")

spp <- species[ , 5:27]
env <- species[ , 1:4]

env <- env %>% unite("VegYr", Year:VegType, remove = FALSE)

env.unin <- species %>% #rename the factors
  mutate(VegType = fct_recode(VegType,
                              "Remnant" = "Emergent",
                              "Remnant" = "Meadow")) 



env.unin <- env.unin %>% unite("VegYr", Year:VegType, remove = FALSE)
str(env.unin)



###### Beta Diversity with two groups #########

env.unin

Year2014.un <- env.unin %>% filter(Year == "2014")
Year2015.un  <- env.unin %>% filter(Year == "2015")

spp2014.un <- Year2014.un[ ,6:28] 
env2014.un <- Year2014.un[ ,1:5]

spp2015.un <- Year2015.un[ ,6:28] # than these sites
env2015.un <- Year2015.un[ ,1:5]


# seperate the groups

uninvaded14 <- Year2014.un %>% filter(VegType == "Remnant")
invaded14 <- Year2014.un %>% filter(VegType == "Invaded")

uninvaded15 <- Year2015.un %>% filter(VegType == "Remnant")
invaded15 <- Year2015.un %>% filter(VegType == "Invaded")


### just species

uninvaded.14 <- uninvaded14[ , 6:28]
uninvaded.15 <- uninvaded15[ , 6:28]

invaded.14 <- invaded14[ ,6:28]
invaded.15<- invaded15[ ,6:28]

## convert to presence absence

uninvaded.14[uninvaded.14 > 0] <- 1
uninvaded.15[uninvaded.15 > 0] <- 1


invaded.14[invaded.14 > 0] <- 1
invaded.15[invaded.15 > 0] <- 1

### Beta Diversity Uninvaded ####

beta.uninvaded <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(4)))
colnames(beta.uninvaded) <- c("Turnover","Nestedness","Sum")


uninvaded14.bd <- beta.multi(uninvaded.14, index.family = "sorensen")
uninvaded14.bd
beta.uninvaded[1,] <- data.frame(matrix(unlist(uninvaded14.bd), nrow = length(1), byrow = T))

uninvaded15.bd <- beta.multi(uninvaded.15, index.family = "sorensen")
uninvaded15.bd
beta.uninvaded[2,] <- data.frame(matrix(unlist(uninvaded15.bd), nrow = length(1), byrow = T))


invaded14.bd <- beta.multi(invaded.14, index.family = "sorensen")
invaded14.bd
beta.uninvaded[3,] <- data.frame(matrix(unlist(invaded14.bd), nrow = length(1), byrow = T))


invaded15.bd <- beta.multi(invaded.15, index.family = "sorensen")
invaded15.bd
beta.uninvaded[4,] <- data.frame(matrix(unlist(invaded15.bd), nrow = length(1), byrow = T))

beta.uninvaded$Vegetation <- as.factor(c("Remnant", "Remnant", "Invaded", "Invaded"))
beta.uninvaded$Year <- as.factor(c("2014","2015","2014","2015"))  
beta.uninvaded$N <- c("14","14","6","6") 


write.csv(beta.uninvaded, "Data/beta_diversity_twogroups.csv")

#### Beta Diversity Figure #####

## Turnover
turn.figu <- ggplot(beta.uninvaded, aes(x = Vegetation, y = Turnover, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Turnover"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  scale_fill_manual(values = c("#fc8d62","#35978f")) 


turn.figu  


# Nestedness

nest.figu <- ggplot(beta.uninvaded, aes(x = Vegetation, y = Nestedness, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Nestedness"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) 


nest.figu 

# sum

sum.figu <- ggplot(beta.uninvaded, aes(x = Vegetation, y = Sum, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Sum"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#fc8d62","#35978f"))

sum.figu

# put all three together

order.un <- ggarrange(nest.figu, turn.figu, sum.figu,
                      ncol = 3, common.legend = TRUE, legend = "right")
order.un

ggsave("Figures/BetaDiversity_uninvaded_invaded.jpeg", order.un)



######### Null Model ############

# combine both years for each veg type for gamma

###### Null with all Species (veg x year) #######

#### Uninvaded  ############

species.uninvaded <- env.unin %>% filter(VegType == "Remnant") # only uninvaded sites
species.unv <- species.uninvaded[ ,6:28] # remove categorical variables
species.unv[species.unv > 0] <- 1 # make pres/abs

unv.year <- species.uninvaded[ , 3] # take out year to add back in

# Create a results data frame for each veg x year data

uninvaded.null.14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(uninvaded.null.14) <- c("Turnover","Nestedness","Sum")

uninvaded.null.15 <-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(uninvaded.null.15) <- c("Turnover","Nestedness","Sum")


# run it 1000 times 

for (i in 1:1000){ #for 1000 iterations
  tempm <- as.data.frame(randomizeMatrix(species.unv, null.model = "independentswap")) #Randomize the data
  tempm[,24] <- unv.year # add year back in 
  un14v <- tempm %>% filter(V24 == "2014") # seperate the matrix based on year
  un15v <- tempm %>% filter(V24 == "2015")
  un14bv <- beta.multi(un14v[,1:23], index.family = "sorensen") # Select just species and calculate beta diversity/nest/turnover
  un15bv <- beta.multi(un15v[,1:23], index.family = "sorensen")
  uninvaded.null.14[i,] <- data.frame(matrix(unlist(un14bv), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  uninvaded.null.15[i,] <- data.frame(matrix(unlist(un15bv), nrow = length(1), byrow = T))
}


## summarize results from null model 
uninvaded.2014 <- uninvaded.null.14 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(14)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(14)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(14)))

uninvaded.2014

# uninvaded 2015
uninvaded.2015 <- uninvaded.null.15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(14)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(14)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(14)))

uninvaded.2015

uninvaded <- rbind(uninvaded.2014, uninvaded.2015)  # combine both years
uninvaded$Year <- as.factor(c("2014","2015")) # add year
uninvaded$Vegetation <- as.factor(c("Remnant", "Remnant")) # add veg
uninvaded


###### Invaded Null model ########

species.invaded <- env.unin %>% filter(VegType == "Invaded") # only uninvaded sites
species.inv <- species.invaded[ ,6:28] # remove categorical variables
species.inv[species.inv > 0] <- 1 # make pres/abs

inv.year <- species.invaded[ , 3] # take out year to add back in

# Create a results data frame for each veg x year data

invaded.null.14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(invaded.null.14) <- c("Turnover","Nestedness","Sum")

invaded.null.15 <-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(invaded.null.15) <- c("Turnover","Nestedness","Sum")


# run it 1000 times 

for (i in 1:1000){ #for 1000 iterations
  tempin <- as.data.frame(randomizeMatrix(species.inv, null.model = "independentswap")) #Randomize the data
  tempin[,24] <- inv.year # add year back in 
  in14v <- tempin %>% filter(V24 == "2014") # seperate the matrix based on year
  in15v <- tempin %>% filter(V24 == "2015")
  in14bv <- beta.multi(in14v[,1:23], index.family = "sorensen") # Select just species and calculate beta diversity/nest/turnover
  in15bv <- beta.multi(in15v[,1:23], index.family = "sorensen")
  invaded.null.14[i,] <- data.frame(matrix(unlist(in14bv), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  invaded.null.15[i,] <- data.frame(matrix(unlist(in15bv), nrow = length(1), byrow = T))
}


## summarize results from null model 
invaded.2014 <- invaded.null.14 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(6)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(6)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(6)))

invaded.2014

# Meadow 2015
invaded.2015 <- invaded.null.15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(6)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(6)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(6)))

invaded.2015

invaded <- rbind(invaded.2014, invaded.2015)  # combine both years
invaded$Year <- as.factor(c("2014","2015")) # add year
invaded$Vegetation <- as.factor(c("Invaded", "Invaded")) # add veg
invaded


########## Put all the data together ###########

unin.null <- rbind(uninvaded, invaded)
unin.null # but together all the results 

unin.null$Turn <- beta.uninvaded$Turnover # add in the real data (from `Diversity measures`)
unin.null$Nest <- beta.uninvaded$Nestedness
unin.null$Sum <- beta.uninvaded$Sum

colnames(unin.null) # reorder for figure
str(unin.null)

unin.null <- unin.null %>% mutate(Vegetation = fct_relevel(Vegetation,
                                              "Invaded", "Remnant"))

write.csv(unin.null, "Data/Beta_null_true_invunin.csv")


###### Figures, points for real data and null model error bars ########


# sum
sum.figure <- ggplot(unin.null, aes(x = Vegetation, y = Sum, colour = Year, shape = Year)) +
  geom_point(position = position_dodge(0.6), size = 5) +
  geom_errorbar(aes(ymin = S.avg - S.CI, ymax = S.avg + S.CI),
                colour = "black",
                size = 1,
                width = 0.3,
                position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Beta Diversity Sum"))) +
  theme_classic() +
  ylim(0, 1) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_colour_manual(values = c("#fc8d62","#35978f"))

sum.figure


# nestedness 
nest.figure <- ggplot(unin.null, aes(x = Vegetation, y = Nest, colour = Year, shape = Year)) +
  geom_point(position = position_dodge(0.6), size = 5) +
  geom_errorbar(aes(ymin = N.avg - N.CI, ymax = N.avg + N.CI),
                colour = "black",
                size = 1,
                width = 0.3,
                position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Nestedness"))) +
  theme_classic() +
  ylim(0, 1) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_colour_manual(values = c("#fc8d62","#35978f"))

nest.figure

## turnover
turn.figure <- ggplot(unin.null, aes(x = Vegetation, y = Turn, colour = Year, shape = Year)) +
  geom_point(position = position_dodge(0.6), size = 5) +
  geom_errorbar(aes(ymin = T.avg - T.CI, ymax = T.avg + T.CI),
                colour = "black",
                size = 1,
                width = 0.3,
                position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Turnover"))) +
  theme_classic() +
  ylim(0, 1) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_colour_manual(values = c("#fc8d62","#35978f"))

turn.figure


null.points.unin <- ggarrange(sum.figure, turn.figure, nest.figure,
                         ncol = 3, common.legend = TRUE, legend = "bottom",
                         labels = "AUTO",
                         hjust = -7,
                         vjust = 2.5)

null.points.unin

ggsave("Figures/beta_null_true_points_unin.TIFF", null.points.unin,
       width = 9.3, height = 5.2, dpi = 150, units = "in")


