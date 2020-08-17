
library(vegan)
library(tidyverse)
library(ggpubr)

### Import Bird matrix ####


species <- read.csv("Data/Species matrix_raw.csv")
species

###### Null Models for Nestedness/Turnover/Beta #########

library(betapart) # beta diversity calculations
library(picante) # randomizing matrix


# combine both years for each veg type for gamma

###### Null with all Species (veg x year) #######

#### Meadow null model  ############

species.med <- species %>% filter(VegType == "Meadow") # only meadow sites
species.m <- species.med[,5:27] # remove categorical variables
species.m[species.m > 0] <- 1 # make pres/abs

m.cat <- species.med[,2] # take out year to add back in

# Create a results data frame for each veg x year data

med14v <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(med14v) <- c("Turnover","Nestedness","Sum")

med15v <-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(med15v) <- c("Turnover","Nestedness","Sum")



for (i in 1:1000){ #for 1000 iterations
  tempm <- as.data.frame(randomizeMatrix(species.m, null.model = "independentswap")) #Randomize the data
  tempm[,24] <- m.cat # add year back in 
  md14v <- tempm %>% filter(V24 == "2014") # seperate the matrix based on year
  md15v <- tempm %>% filter(V24 == "2015")
  med14bv <- beta.multi(md14v[,1:23], index.family = "sorensen") # Select just species and calculate beta diversity/nest/turnover
  med15bv <- beta.multi(md15v[,1:23], index.family = "sorensen")
  med14v[i,] <- data.frame(matrix(unlist(med14bv), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  med15v[i,] <- data.frame(matrix(unlist(med15bv), nrow = length(1), byrow = T))
}

# Meadow 2014

## summarize results from null model 
mveg.2014 <- med14v %>% 
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
mveg.2014

# Meadow 2015
mveg.2015 <- med15v %>% 
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
mveg.2015

meadow <- rbind(mveg.2014, mveg.2015)  # combine both years
meadow$Year <- as.factor(c("2014","2015")) # add year
meadow$Vegetation <- as.factor(c("Meadow", "Meadow")) # add veg
meadow


###### Invaded Null model ########

inv14v <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(inv14v) <- c("Turnover","Nestedness","Sum")

inv15v <-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(inv15v) <- c("Turnover","Nestedness","Sum")

## create data frame

species.inv <- species %>% filter(VegType == "Invaded") # only phrag sites
species.i <- species.inv[,5:27] # remove categorical variables
species.i[species.i > 0] <- 1 # make pres/abs

i.cat <- species.inv[,2] # take out the year
i.cat


for (i in 1:1000){ #for 500 iterations
  tempi <- as.data.frame(randomizeMatrix(species.i, null.model = "independentswap")) #Randomize the data
  tempi[,24] <- i.cat
  in14v <- tempi %>% filter(V24 == "2014")
  in15v <- tempi %>% filter(V24 == "2015")
  in14bv <- beta.multi(in14v[,1:23], index.family = "sorensen") #Run the beta diversity for each year's randomized data
  in15bv <- beta.multi(in15v[,1:23], index.family = "sorensen")
  inv14v[i,] <- data.frame(matrix(unlist(in14bv), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  inv15v[i,] <- data.frame(matrix(unlist(in15bv), nrow = length(1), byrow = T))
}

# Invaded 2014
iveg.2014 <- inv14v %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(6)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(6)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(6))) # sample size is 6
iveg.2014

# Invaded 2015
iveg.2015 <- inv15v %>% 
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
iveg.2015

invaded <- rbind(iveg.2014, iveg.2015)
invaded$Year <- as.factor(c("2014","2015"))
invaded$Vegetation <- as.factor(c("Invaded", "Invaded"))

invaded

###### Emergent Null model ########

emg14v <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(emg14v) <- c("Turnover","Nestedness","Sum")

emg15v <-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(1000)))
colnames(emg15v) <- c("Turnover","Nestedness","Sum")

## create data frame

species.emg <- species %>% filter(VegType == "Emergent") # only emergent sites
species.e <- species.emg[,5:27] # remove categorical variables
species.e[species.e > 0] <- 1 # make pres/abs

e.cat <- species.emg[,2] # take out the year
e.cat


for (i in 1:1000){ #for 500 iterations
  tempe <- as.data.frame(randomizeMatrix(species.e, null.model = "independentswap")) #Randomize the data
  tempe[,24] <- e.cat
  em14v <- tempe %>% filter(V24 == "2014")
  em15v <- tempe %>% filter(V24 == "2015")
  em14bv <- beta.multi(em14v[,1:23], index.family = "sorensen") #Run the beta diversity for each year's randomized data
  em15bv <- beta.multi(em15v[,1:23], index.family = "sorensen")
  emg14v[i,] <- data.frame(matrix(unlist(em14bv), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  emg15v[i,] <- data.frame(matrix(unlist(em15bv), nrow = length(1), byrow = T))
}

# Emergent 2014
eveg.2014 <- emg14v %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))
eveg.2014

# Emergent 2015
eveg.2015 <- emg15v %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))
eveg.2015

emergent <- rbind(eveg.2014, eveg.2015)
emergent$Year <- as.factor(c("2014","2015"))
emergent$Vegetation <- as.factor(c("Emergent", "Emergent"))
emergent

########## Put all the data together ###########
reference <- rbind(emergent, meadow)
all.null <- rbind(reference, invaded)
all.null # but together all the results 

all.null$Turn <- as.numeric(c("0.481","0.446","0.406","0.563","0.395","0.281")) # add in the real data (from `Diversity measures`)
all.null$Nest <- as.numeric(c("0.132","0.143","0.1","0.094","0.183","0.151"))
all.null$Sum <- as.numeric(c("0.614","0.589","0.506","0.656","0.578","0.432"))

all.null$Vegetation <- factor(all.null$Vegetation, levels = c("Meadow","Invaded", "Emergent")) # reorder for figure

write.csv(all.null, "Data/Beta_null_true.csv")


###### Figures, points for real data and null model error bars ########

str(all.null)

# sum
sum.figure <- ggplot(all.null, aes(x = Vegetation, y = Sum, colour = Year, shape = Year)) +
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
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

sum.figure


# nestedness 
nest.figure <- ggplot(all.null, aes(x = Vegetation, y = Nest, colour = Year, shape = Year)) +
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
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

nest.figure

## turnover
turn.figure <- ggplot(all.null, aes(x = Vegetation, y = Turn, colour = Year, shape = Year)) +
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
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

turn.figure


null.points <- ggarrange(sum.figure, nest.figure, turn.figure,
                         ncol = 3, common.legend = TRUE, legend = "bottom")

null.points

ggsave("Figures/beta_null_true_points.jpg", null.points)



######### testing with gamma as the whole dataset
# Don't need to do this because there is an interaction but tried it anyway

###### Null with all Species (veg x year) #######
species.pa <- species[,5:27]
species.pa[species.pa > 0] <- 1


# Create a results data frame for each veg x year data

med14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(med14) <- c("Turnover","Nestedness","Sum")

med15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(med15) <- c("Turnover","Nestedness","Sum")

inv14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(inv14) <- c("Turnover","Nestedness","Sum")

inv15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(inv15) <- c("Turnover","Nestedness","Sum")

emg14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(emg14) <- c("Turnover","Nestedness","Sum")

emg15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(emg15) <- c("Turnover","Nestedness","Sum")

# loop beta diversity through each matrix 

temp.cat <- species[,2:3]
test <- species.pa 
  
test[,24:26] <- temp.cat



for (i in 1:500){ #for 500 iterations
  temp <- as.data.frame(randomizeMatrix(species.pa, null.model = "independentswap")) #Randomize the data
  temp[,24:25] <- temp.cat
  md14 <- temp %>% filter(VegType == "Meadow" & Year == "2014")
  md15 <- temp %>% filter(VegType == "Meadow" & Year == "2015")
  in14 <- temp %>% filter(VegType == "Invaded" & Year == "2014")
  in15 <- temp %>% filter(VegType == "Invaded" & Year == "2015")
  em14 <- temp %>% filter(VegType == "Emergent" & Year == "2014")
  em15 <- temp %>% filter(VegType == "Emergent" & Year == "2015")
  med14b <- beta.multi(md14[,1:23], index.family = "sorensen") #Run the beta diversity for each year's randomized data
  med15b <- beta.multi(md15[,1:23], index.family = "sorensen")
  inv14b <- beta.multi(in14[,1:23], index.family = "sorensen") #Run the beta diversity for each year's randomized data
  inv15b <- beta.multi(in15[,1:23], index.family = "sorensen")
  emg14b <- beta.multi(em14[,1:23], index.family = "sorensen") #Run the beta diversity for each year's randomized data
  emg15b <- beta.multi(em15[,1:23], index.family = "sorensen")
  med14[i,] <- data.frame(matrix(unlist(med14b), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  med15[i,] <- data.frame(matrix(unlist(med15b), nrow = length(1), byrow = T))
  inv14[i,] <- data.frame(matrix(unlist(inv14b), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  inv15[i,] <- data.frame(matrix(unlist(inv15b), nrow = length(1), byrow = T))
  emg14[i,] <- data.frame(matrix(unlist(emg14b), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  emg15[i,] <- data.frame(matrix(unlist(emg15b), nrow = length(1), byrow = T))
}

# Meadow 2014
m.2014 <- med14 %>% 
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
m.2014

# Meadow 2015
m.2015 <- med15 %>% 
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
m.2015

meadow <- rbind(m.2014, m.2015)

# Invaded 2014
i.2014 <- inv14 %>% 
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
i.2014

# Invaded 2015
i.2015 <- inv15 %>% 
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
i.2015

invaded <- rbind(e.2014, e.2015)

# Emergent 2014
e.2014 <- emg14 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))
e.2014

# Emergent 2014
e.2015 <- emg15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))
e.2015

emergent <- rbind(e.2014, e.2015)

reference <- rbind(emergent, meadow)
all.null <- rbind(reference, invaded)

all.null$Vegetation <- as.factor(c("Emergent","Emergent","Meadow","Meadow","Invaded","Invaded"))
all.null$Year <- as.factor(c("2014","2015","2014","2015","2014","2015"))
all.null$Turn <- as.numeric(c("0.481","0.446","0.406","0.563","0.395","0.281"))
all.null$Nest <- as.numeric(c("0.132","0.143","0.1","0.094","0.183","0.151"))
all.null$Sum <- as.numeric(c("0.614","0.589","0.506","0.656","0.578","0.432"))
all.null$Vegetation <- factor(all.null$Vegetation, levels = c("Meadow","Invaded", "Emergent")) 

colnames(all.null)
str(all.null)

all.null

## Turnover

turn.fig <- ggplot(all.null, aes(x = Vegetation, y = Turn, colour = Year, shape = Year)) +
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
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

turn.fig

### Sum 

sum.fig <- ggplot(all.null, aes(x = Vegetation, y = Sum, colour = Year, shape = Year)) +
  geom_point(position = position_dodge(0.6), size = 5) +
  geom_errorbar(aes(ymin = S.avg - S.CI, ymax = S.avg + S.CI),
                colour = "black",
                size = 1,
                width = 0.3,
                position = position_dodge(0.6)) +
  labs(x = " ",
       y = expression(paste("Sum"))) +
  theme_classic() +
  ylim(0, 1) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

sum.fig


# nestedness 
nest.fig <- ggplot(all.null, aes(x = Vegetation, y = Nest, colour = Year, shape = Year)) +
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
  scale_colour_manual(values = c("#fc8d62","#1f78b4"))

nest.fig



null.fig.panel <- ggarrange(sum.fig, nest.fig, turn.fig,
                         ncol = 3, common.legend = TRUE, legend = "bottom")

ggarrange(null.fig.panel, null.points,
          nrow = 2,
          labels = "AUTO")