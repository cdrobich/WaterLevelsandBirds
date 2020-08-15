
library(vegan)
library(tidyverse)
library(ggpubr)

### Import Bird matrix ####

env <- read.csv("Data/Env_variables.csv")
species <- read.csv("Data/Species matrix_raw.csv")

###### Null Models for Nestedness/Turnover/Beta #########
library(betapart)
library(picante)

### someone else's code that is pretty

do.call(rbind, lapply(1:2, function(i) {
  x <-  randomizeMatrix(com.dat, null.model = "richness", iterations = 1000)
  unlist(beta.multi.abund(x, index.family = "bray"))
}))

# a random value & confidence intervals - 
# beta of each of 500 examples & take differences & CI 

# combine both years bc otherwise not random enough
# Sites 1 - 6 are 2014 and 7 - 12 are 2015

###### Meadow Null #######
Meadow <- species %>% filter(VegType== "Meadow")
Meadow <- Meadow[ ,5:27]
Meadow[Meadow > 0] <- 1


# Create a results data frame for each year's data

bd14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(bd14) <- c("Turnover","Nestedness","Sum")

bd15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(bd15) <- c("Turnover","Nestedness","Sum")

# loop beta diversity through each matrix 

for (i in 1:500){ #for 500 iterations
  temp <- randomizeMatrix(Meadow, null.model = "independentswap") #Randomize the data
  veg14 <- temp[1:6,] #split the data into years
  veg15 <- temp[7:12,]
  b14 <- beta.multi(veg14, index.family = "sorensen") #Run the beta diversity for each year's randomized data
  b15 <- beta.multi(veg15, index.family = "sorensen")
  bd14[i,] <- data.frame(matrix(unlist(b14), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  bd15[i,] <- data.frame(matrix(unlist(b15), nrow = length(1), byrow = T))
}

# Meadow 2014
sum.2014 <- bd14 %>% 
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
sum.2014

# Meadow 2015 
sum.2015 <- bd15 %>% 
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
sum.2015

meadow45 <- rbind(sum.2014, sum.2015)
meadow45$Veg <- c("Meadow","Meadow")
meadow45$Year <- c("2014","2015")
meadow45

##### Invaded Null models #########

Invaded <- species %>% filter(VegType== "Invaded")
Invaded <- Invaded[ ,5:27]
Invaded[Invaded > 0] <- 1

# data frames
inv14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(inv14) <- c("Turnover","Nestedness","Sum")

inv15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(inv15) <- c("Turnover","Nestedness","Sum")

# loop

for (i in 1:500){ #for 500 iterations
  temp <- randomizeMatrix(Invaded, null.model = "independentswap") #Randomize the data
  veg14 <- temp[1:6,] #split the data into years
  veg15 <- temp[7:12,]
  b14 <- beta.multi(veg14, index.family = "sorensen") #Run the beta diversity for each year's randomized data
  b15 <- beta.multi(veg15, index.family = "sorensen")
  inv14[i,] <- data.frame(matrix(unlist(b14), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  inv15[i,] <- data.frame(matrix(unlist(b15), nrow = length(1), byrow = T))
}


# Invaded 2014
inv.2014 <- inv14 %>% 
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
inv.2014

# Invaded 2015 
inv.2015 <-inv15 %>% 
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
inv.2015

inv45 <- rbind(inv.2014, inv.2015)
inv45$Veg <- c("Invaded","Invaded")
inv45$Year <- c("2014","2015")
inv45

##### Emergent Null models #########

Emerg <- species %>% filter(VegType== "Emergent")
Emerg <- Emerg[ ,5:27]
Emerg[Emerg > 0] <- 1

# data frames
emg14 <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(emg14) <- c("Turnover","Nestedness","Sum")

emg15<-data.frame(matrix(as.numeric(0), ncol=(3), nrow=(500)))
colnames(emg15) <- c("Turnover","Nestedness","Sum")

# loop

for (i in 1:500){ #for 500 iterations
  temp <- randomizeMatrix(Emerg, null.model = "independentswap") #Randomize the data
  veg14 <- temp[1:6,] #split the data into years
  veg15 <- temp[7:12,]
  b14 <- beta.multi(veg14, index.family = "sorensen") #Run the beta diversity for each year's randomized data
  b15 <- beta.multi(veg15, index.family = "sorensen")
  emg14[i,] <- data.frame(matrix(unlist(b14), nrow = length(1), byrow = T)) #Put the beta results in the correct results df
  emg15[i,] <- data.frame(matrix(unlist(b15), nrow = length(1), byrow = T))
}


# Emergent 2014
emg.2014 <- emg14 %>% 
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
emg.2014

# Emergent 2015 
emg.2015 <- emg15 %>% 
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
emg.2015

emerg45 <- rbind(emg.2014, emg.2015)
emerg45$Veg <- c("Emergent","Emergent")
emerg45$Year <- c("2014","2015")
emerg45


Beta.null <- full_join(meadow45, inv45)
Beta.null <- full_join(Beta.null, emerg45)

Beta.null.sum <- Beta.null %>% 
  mutate(N.up = (N.avg + N.CI),
         N.low = (N.avg - N.CI),
         T.up = (T.avg + T.CI),
         T.low = (T.avg - T.CI),
         S.up = (S.avg + S.CI),
         S.low = (S.avg - S.CI))

Beta.null.sum

write.csv(Beta.null.sum, "Data/beta_null_summary.csv")

Beta.null.sum$Year <- as.factor(Beta.null.sum$Year)
Beta.null.sum$Veg <- factor(Beta.null.sum$Veg , levels = c("Meadow","Invaded", "Emergent")) 

# Nestedness null 
nest.null <- ggplot(Beta.null.sum, aes(y = N.avg, x = Veg, fill = Year)) +
  geom_col(color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = N.avg - N.CI, ymax = N.avg + N.CI),
                width = 0.2,
                position = position_dodge(0.9)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#bdbdbd","#636363")) +
  labs(x = " ",
       y = expression(paste("Nestedness (null)"))) +
  ylim(0, 0.8)

nest.null

# Turnover null
turn.null <- ggplot(Beta.null.sum, aes(y = T.avg, x = Veg, fill = Year)) +
  geom_col(color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = T.avg - T.CI, ymax = T.avg + T.CI),
                width = 0.2,
                position = position_dodge(0.9)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#bdbdbd","#636363")) +
  labs(x = " ",
       y = expression(paste("Turnover (null)"))) +
  ylim(0, 0.8)

turn.null



# Sum null
sum.null <- ggplot(Beta.null.sum, aes(y = S.avg, x = Veg, fill = Year)) +
  geom_col(color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = S.avg - S.CI, ymax = S.avg + S.CI),
                width = 0.2,
                position = position_dodge(0.9)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#bdbdbd","#636363")) +
  labs(x = " ",
       y = expression(paste("Sum (null)"))) +
  ylim(0, 0.8)

sum.null


order.null <- ggarrange(nest.null, turn.null, sum.null,
                        ncol = 3, common.legend = TRUE, legend = "right")
order.null

panel <- ggarrange(order, order.null,
                   nrow = 2)
panel

ggsave("Figures/Beta_null_true.jpg", panel)


###### Points and Null CI


true.null <- read.csv("Data/betadiversity_true_null.csv")
true.null
true.null$Year <- as.factor(true.null$Year)
true.null$Vegetation <- factor(true.null$Vegetation , levels = c("Meadow","Invaded", "Emergent")) 

# sum
sum.figure <- ggplot(true.null, aes(x = Vegetation, y = Sum, colour = Year, shape = Year)) +
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

sum.figure


# nestedness 
nest.figure <- ggplot(true.null, aes(x = Vegetation, y = Nest, colour = Year, shape = Year)) +
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
turn.figure <- ggplot(true.null, aes(x = Vegetation, y = Turn, colour = Year, shape = Year)) +
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

all.null$Vegetation <- c("Emergent","Emergent","Meadow","Meadow","Invaded","Invaded")
all.null$Year <- as.factor(c("2014","2015","2014","2015","2014","2015"))
all.null$Turn <- as.numeric(c("0.481","0.446","0.406","0.563","0.395","0.281"))
all.null$Nest <- as.numeric(c("0.132","0.143","0.1","0.094","0.183","0.151"))
all.null$Sum <- as.numeric(c("0.614","0.589","0.506","0.656","0.578","0.432"))

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