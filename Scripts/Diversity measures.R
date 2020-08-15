
library(vegan)
library(tidyverse)
library(ggpubr)

### Import Bird matrix ####

env <- read.csv("Data/Env_variables.csv")
species <- read.csv("Data/Species matrix_raw.csv")

spp <- species[ , 5:27]

beta <- betadiver(spp, "-1")
mod <- with(env, betadisper(beta, VegYr))

anova(mod)

boxplot(mod,
        xlab = " ",
        ylab = "Distance to Centroid (Beta Diversity)")
       


#### Alpha diversity ####
species$Year <- as.factor(species$Year)
str(species)

H <- diversity(spp) # Shannon index
species$Shannon <- H 

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


####### Beta Diversity #######

library(adespatial)
library(betapart)
citation("betapart")

Year2014 <- species %>% filter(Year == "2014")
Year2015 <- species %>% filter(Year == "2015")

spp2014 <- Year2014[ ,5:27]
env2014 <- Year2014[ ,1:4]

spp2015 <- Year2015[ ,5:27]
env2015 <- Year2015[ ,1:4]


#### Create Bray-Curtis dissimilarity matrix ####

spp2014bc <- vegdist(spp2014, method = "bray",
                     binary = FALSE)
spp2014bc


spp2015bc <- vegdist(spp2015, method = "bray",
                     binary = FALSE)
spp2015bc

## beta diversity 

## beta diversity and nestedness and turnover by veg x year 

# seperate the groups

meadow14 <- Year2014 %>% filter(VegType == "Meadow")
phrag14 <- Year2014 %>% filter(VegType == "Invaded")
emerg14 <- Year2014 %>% filter(VegType == "Emergent")

meadow15 <- Year2015 %>% filter(VegType == "Meadow")
phrag15 <- Year2015 %>% filter(VegType == "Invaded")
emerg15 <- Year2015 %>% filter(VegType == "Emergent")

### just species

spp2014 <- Year2014[ ,5:27]

meadow14sp <- meadow14[ , 5:27]
meadow15sp <- meadow15[ , 5:27]

phrag14sp <- phrag14[ ,5:27]
phrag15sp <- phrag15[ ,5:27]

emerg14sp <- emerg14[ ,5:27]
emerg15sp <- emerg15[ ,5:27]

## convert to presence absence

meadow14sp[meadow14sp > 0] <- 1
meadow15sp[meadow15sp > 0] <- 1

phrag14sp[phrag14sp > 0] <- 1
phrag15sp[phrag15sp > 0] <- 1

emerg14sp[emerg14sp > 0] <- 1
emerg15sp[emerg15sp > 0] <- 1


### Nestedness and Turnover ####

beta.multi(meadow14sp, index.family = "sorensen")
#$beta.SIM - Turnover
#[1] 0.40625

#$beta.SNE - Nestedness
#[1] 0.1002435

#$beta.SOR - Overall value
#[1] 0.5064935

beta.multi(meadow15sp, index.family = "sorensen")

#$beta.SIM - Turnover
#[1] 0.5625

#$beta.SNE - Nestedness
#[1] 0.09394172

#$beta.SOR - Sum
#[1] 0.6564417

beta.multi(phrag14sp, index.family = "sorensen")

#$beta.SIM - Turnover
#[1] 0.3947368

#$beta.SNE - Nestedness
#[1] 0.1832448

#$beta.SOR - Sum
#[1] 0.5779817

beta.multi(phrag15sp, index.family = "sorensen")

#$beta.SIM
#[1] 0.28125

#$beta.SNE
#[1] 0.1508488

#$beta.SOR
#[1] 0.4320988

beta.multi(emerg14sp, index.family = "sorensen")
#$beta.SIM
#[1] 0.4814815

#$beta.SNE
#[1] 0.1323116

#$beta.SOR
#[1] 0.6137931

beta.multi(emerg15sp, index.family = "sorensen")

#$beta.SIM
#[1] 0.4464286

#$beta.SNE
#[1] 0.1429754

#$beta.SOR
#[1] 0.589404


## made an excel file bc I am lazy; saved in data/betadiversity_output
### turnover nestedness figures ####

beta <- read.csv("Data/betadiversity_output.csv")
beta$Year <- as.factor(beta$Year)
beta$Vegetation <- factor(beta$Vegetation, levels = c("Meadow","Invaded", "Emergent")) 

## Turnover
turn.fig <- ggplot(beta, aes(x = Vegetation, y = Turn, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Turnover"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) 


turn.fig  


# Nestedness

nest.fig <- ggplot(beta, aes(x = Vegetation, y = Nest, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Nestedness"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#fc8d62","#1f78b4")) 


nest.fig 

# sum

sum.fig <- ggplot(beta, aes(x = Vegetation, y = Sum, fill = Year)) +
  geom_col(color = "black", position = position_dodge()) +
  labs(x = " ",
       y = expression(paste("Sum"))) +
  theme_classic() +
  ylim(0, 0.8) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = c("#fc8d62","#1f78b4"))

sum.fig


order <- ggarrange(nest.fig, turn.fig, sum.fig,
                   ncol = 3, common.legend = TRUE, legend = "right")
order

ggsave("Figures/BetaDiversity_partition.jpeg", order)


###### Null Models for Nestedness/Turnover/Beta #########

library(picante)

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
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
sum.2014

# Meadow 2015 
sum.2015 <- bd15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
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
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
inv.2014

# Invaded 2015 
inv.2015 <-inv15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
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
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
emg.2014

# Emergent 2015 
emg.2015 <- emg15 %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(n)),
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(n)),
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(n)))
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
       y = expression(paste("Sum Beta Diversity (null)"))) +
  ylim(0, 0.8)

sum.null


order.null <- ggarrange(nest.null, turn.null, sum.null,
                   ncol = 3, common.legend = TRUE, legend = "right")
order.null

panel <- ggarrange(order, order.null,
          nrow = 2)

ggsave("Figures/Beta_null_true.jpg", panel)

########### Local Contribution to Beta Diversity ############

BD2014 <- beta.div(spp2014bc, method = "percentdiff",
         sqrt.D = FALSE, samp = FALSE,
         nperm = 999)

BD2014

Year2014$LCBD <- BD2014$LCBD

BD2014$LCBD[BD2014$LCBD > mean(BD2014$LCBD)] # LCBD > average

#4           5           6           7           9          10         13         19 
#0.12124900  0.06212209  0.05502379  0.07485454  0.05221660 0.05759652 0.05109404 0.05577960 

YearF <- ggplot(Year2014, aes(y = LCBD, x = VegType)) +
  geom_point() +
  theme_classic()

YearF

#Beta
# max value of BD is 0.5

#SStotal      BDtotal 
#0.44326161   0.02332956

#$LCBD; each sites local contribution to BD
#1          2          3          4 
#0.04645560 0.04315607 0.02153647 0.12124900 
#5          6          7          8 
#0.06212209 0.05502379 0.07485454 0.04894537 
#9         10         11         12 
#0.05221660 0.05759652 0.04786658 0.03867879 
#13         14         15         16 
#0.05109404 0.02905092 0.04796502 0.04626943 
#17         18         19         20 
#0.02844023 0.04150518 0.05577960 0.03019418 


BD2015 <- beta.div(spp2015bc, method = "percentdiff",
         sqrt.D = FALSE, samp = FALSE,
         nperm = 999)

Year2015$LCBD <- BD2015$LCBD

BD2015$LCBD[BD2015$LCBD > mean(BD2015$LCBD)] #which sites have a LCBD > mean

#1           2           3           7           11          13          16          19          20 
#0.05014927  0.08997777  0.05113195  0.05064635  0.05406081  0.09452751  0.05653058  0.06947543  0.06421346 


Year5 <- ggplot(Year2015, aes(y = LCBD, x = VegType)) +
  geom_point() +
  theme_classic()

# $beta

#SStotal     BDtotal 
#0.41512107 0.02184848 

#LCBD
#1          2          3          4          5          6          7          8          9         10         11 
#0.05014927 0.08997777 0.05113195 0.04146118 0.02640905 0.02824501 0.05064635 0.02984169 0.04263854 0.04429739 0.05406081 
#12         13         14         15         16         17         18         19         20 
#0.04266220 0.09452751 0.03850845 0.03558245 0.05653058 0.04542190 0.04421900 0.06947543 0.06421346 


# Combine the data sets and calculate differences
library(car)
library(agricolae)


both_years <- full_join(Year2014, Year2015)

check <- lm(LCBD ~ VegType*Year, data = both_years)
summary(check)

Anova(check, type = "2")

plot(residuals(check)~fitted(check))


# Type III SS
#VegType      0.0017365  2  2.3904   0.1068    
#Year         0.0000968  1  0.2665   0.6090    
#VegType:Year 0.0001722  2  0.2371   0.7902    
#Residuals    0.0123492 34  

# Type II SS
#Sum Sq Df F value  Pr(>F)  
#VegType      0.0020988  2  2.8893 0.06936 .
#Year         0.0000000  1  0.0000 1.00000  
#VegType:Year 0.0001722  2  0.2371 0.79024  
#Residuals    0.0123492 34 

sum <- both_years %>% group_by(VegType, Year) %>% 
  summarise(Avg.LCBD = mean(LCBD),
            sd.LCBD = sd(LCBD),
            min.LCBD = min(LCBD),
            max.LCBD = max(LCBD))

#VegType  Year  Avg.LCBD sd.LCBD min.LCBD max.LCBD
#<chr>    <fct>    <dbl>   <dbl>    <dbl>    <dbl>
#1 Emergent 2014  0.0610  0.0252    0.0415   0.121 
#2 Emergent 2015  0.0561  0.0197    0.0356   0.0945
#3 Invaded  2014  0.0394  0.0125    0.0215   0.0550
#4 Invaded  2015  0.0440  0.00836   0.0298   0.0541
#5 Meadow   2014  0.0460  0.0167    0.0284   0.0749
#6 Meadow   2015  0.0479  0.0231    0.0264   0.0900

sum2 <- both_years %>% group_by(VegType) %>% 
  summarise(Avg.LCBD = mean(LCBD),
            sd.LCBD = sd(LCBD),
            min.LCBD = min(LCBD),
            max.LCBD = max(LCBD))
sum2

#VegType  Avg.LCBD sd.LCBD  min.LCBD  max.LCBD
#<chr>       <dbl>   <dbl>    <dbl>    <dbl>
#1 Emergent   0.0585  0.0220   0.0356   0.121 
#2 Invaded    0.0417  0.0104   0.0215   0.0550
#3 Meadow     0.0469  0.0193   0.0264   0.0900

library(Hmisc)
colnames(both_years)

LCBD.fig <- ggplot(both_years, aes(x = VegType, y = LCBD)) +
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
       y = expression(paste("LCBD"))) + 
  scale_color_manual(values = c("#fc8d62","#1f78b4")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank")

LCBD.fig

compare <- data.frame(Year2014$Site, Year2014$VegType)

compare$LCBD2014 <- Year2014$LCBD 
compare$LCBD2015 <- Year2015$LCBD 
colnames(compare)

compare2 <- compare %>% mutate(LCBDdiff = LCBD2015 - LCBD2014) 

#Year2014.Site Year2014.VegType   LCBD2014   LCBD2015         diff
#1        CM10_14           Meadow 0.04645560 0.05014927  0.003693677
#2        CM19_14          Invaded 0.04315607 0.08997777  0.046821700
#3         CM2_14          Invaded 0.02153647 0.05113195  0.029595483
#4        CM4R_14         Emergent 0.12124900 0.04146118 -0.079787820
#5        CM5R_14         Emergent 0.06212209 0.02640905 -0.035713042
#6         CM6_14          Invaded 0.05502379 0.02824501 -0.026778776
#7         CM9_14           Meadow 0.07485454 0.05064635 -0.024208188
#8         LP1_14          Invaded 0.04894537 0.02984169 -0.019103683
#9        LP10_14         Emergent 0.05221660 0.04263854 -0.009578064
#10      LP10R_14         Emergent 0.05759652 0.04429739 -0.013299131
#11       LP12_14           Meadow 0.04786658 0.05406081  0.006194228
#12      LP12R_14          Invaded 0.03867879 0.04266220  0.003983412
#13       LP15_14         Emergent 0.05109404 0.09452751  0.043433467
#14      LP16R_14          Invaded 0.02905092 0.03850845  0.009457529
#15        LP5_14           Meadow 0.04796502 0.03558245 -0.012382564
#16       LP5R_14         Emergent 0.04626943 0.05653058  0.010261152
#17        LP6_14           Meadow 0.02844023 0.04542190  0.016981671
#18       LP6R_14         Emergent 0.04150518 0.04421900  0.002713829
#19       LP7R_14         Emergent 0.05577960 0.06947543  0.013695836
#20       LP8R_14           Meadow 0.03019418 0.06421346  0.034019284




#### Temporal Beta Diversity Index ####

TBI1 <- TBI(spp2014, spp2015, method = "%difference", nperm = 999, test.t.perm = FALSE)
TBI1

TBI1$BCD.mat

compare2$TBI <- TBI1$BCD.mat
compare2$p.TBI <- TBI1$p.TBI
compare2

write.csv(compare2, "Data/LCBD_TBI_data.csv")

# C = gain, B = loss; B > C site has lost species between time 1 and 2 (-)


#B/(2A+B+C) C/(2A+B+C) D=(B+C)/(2A+B+C) Change
#  Site.1  0.05714286  0.6000000        0.6571429    +  
#  Site.2  0.32000000  0.3600000        0.6800000    +  
#  Site.3  0.05000000  0.3000000        0.3500000    +  
#  Site.4  0.30555556  0.5277778        0.8333333    +  
#  Site.5  0.14705882  0.3235294        0.4705882    +  
#  Site.6  0.09090909  0.5000000        0.5909091    +  
#  Site.7  0.04761905  0.6666667        0.7142857    +  
#  Site.8  0.32352941  0.1470588        0.4705882    -  
#  Site.9  0.25000000  0.1875000        0.4375000    -  
#  Site.10 0.30000000  0.3000000        0.6000000    0  
#  Site.11 0.14285714  0.1904762        0.3333333    +  
#  Site.12 0.13333333  0.2666667        0.4000000    +  
#  Site.13 0.29032258  0.2580645        0.5483871    -  
#  Site.14 0.27906977  0.1627907        0.4418605    -  
#  Site.15 0.05405405  0.4594595        0.5135135    +  
#  Site.16 0.12121212  0.4545455        0.5757576    +  
#  Site.17 0.12121212  0.3939394        0.5151515    +  
#  Site.18 0.06666667  0.3555556        0.4222222    +  
#  Site.19 0.12244898  0.1428571        0.2653061    +  
#  Site.20 0.15909091  0.6136364        0.7727273    +  

# test that there is no sig. difference between times

#p.TBI
# 0.179 0.153 0.916 0.017 0.641
# 0.303 0.101 0.650 0.730 0.314 
#0.950 0.837 0.439 0.692 0.537 
# 0.371 0.531 0.739 0.984 0.049

# only site 4 and site 20 are significant

## plot TBI
str(compare2)

site.names = rownames(TBI1$BCD.mat)

plot(TBI1, type = "BC", s.names = site.names, 
     pch.loss = 21, pch.gain = 22,
     cex.symb = 4)




