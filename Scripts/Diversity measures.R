
library(vegan)
library(tidyverse)
library(ggpubr)

### Import Bird matrix ####

env <- read.csv("Data/Env_variables.csv")
speciesraw <- read.csv("Data/Species matrix_raw.csv")

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

species <- read.csv("Data/Species matrix_column relativized.csv")

species <- full_join(env, species) # join by site


Year2014 <- species %>% filter(Year == "2014")
Year2015 <- species %>% filter(Year == "2015")

spp2014 <- Year2014[ ,6:28]
env2014 <- Year2014[ ,1:5]

spp2015 <- Year2015[ ,6:28]
env2015 <- Year2015[ ,1:5]


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

spp2014 <- Year2014[ ,6:28]

meadow14sp <- meadow14[ , 6:28]
meadow15sp <- meadow15[ , 6:28]

phrag14sp <- phrag14[ ,6:28]
phrag15sp <- phrag15[ ,6:28]

emerg14sp <- emerg14[ ,6:28]
emerg15sp <- emerg15[ ,6:28]

## convert to presence absence

meadow14sp[meadow14sp > 0] <- 1
meadow15sp[meadow15sp > 0] <- 1

phrag14sp[phrag14sp > 0] <- 1
phrag15sp[phrag15sp > 0] <- 1

emerg14sp[emerg14sp > 0] <- 1
emerg15sp[emerg15sp > 0] <- 1


### Nestedness and Turnover ####

beta <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(6)))
colnames(beta) <- c("Turnover","Nestedness","Sum")


m14.bd <- beta.multi(meadow14sp, index.family = "sorensen")
m14.bd
beta[1,] <- data.frame(matrix(unlist(m14.bd), nrow = length(1), byrow = T))

#$beta.SIM - Turnover
#[1] 0.40625

#$beta.SNE - Nestedness
#[1] 0.1002435

#$beta.SOR - Overall value
#[1] 0.5064935

m15.bd <- beta.multi(meadow15sp, index.family = "sorensen")
m15.bd
beta[2,] <- data.frame(matrix(unlist(m15.bd), nrow = length(1), byrow = T))

#$beta.SIM - Turnover
#[1] 0.5625

#$beta.SNE - Nestedness
#[1] 0.09394172

#$beta.SOR - Sum
#[1] 0.6564417

inv14.bd <- beta.multi(phrag14sp, index.family = "sorensen")
inv14.bd
beta[3,] <- data.frame(matrix(unlist(inv14.bd), nrow = length(1), byrow = T))

#$beta.SIM - Turnover
#[1] 0.3947368

#$beta.SNE - Nestedness
#[1] 0.1832448

#$beta.SOR - Sum
#[1] 0.5779817

inv15.bd <- beta.multi(phrag15sp, index.family = "sorensen")
inv15.bd
beta[4,] <- data.frame(matrix(unlist(inv15.bd), nrow = length(1), byrow = T))


#$beta.SIM
#[1] 0.28125

#$beta.SNE
#[1] 0.1508488

#$beta.SOR
#[1] 0.4320988

emg14.bd <- beta.multi(emerg14sp, index.family = "sorensen")
beta[5,] <- data.frame(matrix(unlist(emg14.bd), nrow = length(1), byrow = T))

#$beta.SIM
#[1] 0.4814815

#$beta.SNE
#[1] 0.1323116

#$beta.SOR
#[1] 0.6137931

emg15.bd <- beta.multi(emerg15sp, index.family = "sorensen")
beta[6,] <- data.frame(matrix(unlist(emg15.bd), nrow = length(1), byrow = T))


#$beta.SIM
#[1] 0.4464286

#$beta.SNE
#[1] 0.1429754

#$beta.SOR
#[1] 0.589404

beta$Vegetation <- as.factor(c("Meadow", "Meadow", "Invaded", "Invaded", "Emergent", "Emergent"))
beta$Year <- as.factor(c("2014","2015","2014","2015","2014","2015"))  
beta$N <- c("6","6","6","6","8","8")  

## Turnover
turn.fig <- ggplot(beta, aes(x = Vegetation, y = Turnover, fill = Year)) +
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

nest.fig <- ggplot(beta, aes(x = Vegetation, y = Nestedness, fill = Year)) +
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

# put all three together

order <- ggarrange(nest.fig, turn.fig, sum.fig,
                   ncol = 3, common.legend = TRUE, legend = "right")
order

ggsave("Figures/BetaDiversity_partition.jpeg", order)


######### Local Contribution to Beta Diversity ############

BD2014 <- beta.div(spp2014bc, method = "percentdiff",
         sqrt.D = FALSE, samp = FALSE,
         nperm = 999)

BD2014

Year2014$LCBD <- BD2014$LCBD

BD2014$LCBD[BD2014$LCBD > mean(BD2014$LCBD)] # LCBD > average

#1           3           4          5          8          9         11         13         15         17         19 
#0.05084490 0.06910410  0.06052732 0.05717447 0.09337527 0.06867763 0.05044497 0.05014559 0.06827220 0.06267009 0.05451351  

YearF <- ggplot(Year2014, aes(y = LCBD, x = VegType)) +
  geom_point() +
  theme_classic()

YearF

#Beta
# max value of BD is 0.5

#SStotal    BDtotal 
#0.25195771 0.01326093 

#$LCBD; each sites local contribution to BD

#1          2          3          4          5          6          7          8          9         10         11         12 
#0.05084490 0.04215959 0.06910410 0.06052732 0.05717447 0.03189931 0.03185241 0.09337527 0.06867763 0.03877596 0.05044497 0.03203308 

#13         14         15         16         17         18         19         20 
#0.05014559 0.04293677 0.06827220 0.02672509 0.06267009 0.03406406 0.05451351 0.03380367


BD2015 <- beta.div(spp2015bc, method = "percentdiff",
         sqrt.D = FALSE, samp = FALSE,
         nperm = 999)

Year2015$LCBD <- BD2015$LCBD

BD2015$LCBD[BD2015$LCBD > mean(BD2015$LCBD)] #which sites have a LCBD > mean

#1          2          5         13         15         17         18         19 
#0.07723281 0.08751646 0.05694796 0.06875478 0.05215225 0.06563237 0.06284933 0.07279255 


Year5 <- ggplot(Year2015, aes(y = LCBD, x = VegType)) +
  geom_point() +
  theme_classic()

Year5

BD2015
# $beta

#SStotal   BDtotal 
#0.2747038 0.0144581 

#1          2          3          4          5          6          7          8          9         10         11         12 
#0.07723281 0.08751646 0.04639343 0.04583418 0.05694796 0.04614452 0.04659774 0.03775481 0.03572251 0.03998418 0.04421550 0.03553016 
#13         14         15         16         17         18         19         20 
#0.06875478 0.01982684 0.05215225 0.03923735 0.06563237 0.06284933 0.07279255 0.01888027 


# Combine the data sets and calculate differences
library(car)
library(agricolae)


both_years <- full_join(Year2014, Year2015)

check <- lm(LCBD ~ VegType*Year, data = both_years)
summary(check)

Anova(check, type = "3")

#Anova Table (Type III tests)

#Response: LCBD
#Sum Sq Df F value Pr(>F)
#(Intercept)  0.0000054  1  0.0175 0.8954
#VegType      0.0007560  2  1.2194 0.3080
#Year         0.0000057  1  0.0183 0.8932
#VegType:Year 0.0007563  2  1.2198 0.3079
#Residuals    0.0105404 34 

plot(residuals(check)~fitted(check))


sum <- both_years %>% group_by(VegType, Year) %>% 
  summarise(Avg.LCBD = mean(LCBD),
            sd.LCBD = sd(LCBD),
            min.LCBD = min(LCBD),
            max.LCBD = max(LCBD))
sum

#VegType   Year Avg.LCBD sd.LCBD min.LCBD max.LCBD
#1 Emergent  2014   0.0488 0.0144    0.0267   0.0687
#2 Emergent  2015   0.0500 0.0216    0.0189   0.0728
#3 Invaded   2014   0.0519 0.0244    0.0319   0.0934
#4 Invaded   2015   0.0400 0.00457   0.0355   0.0466
#5 Meadow    2014   0.0496 0.0147    0.0319   0.0683
#6 Meadow    2015   0.0600 0.0181    0.0458   0.0875

sum2 <- both_years %>% group_by(VegType) %>% 
  summarise(Avg.LCBD = mean(LCBD),
            sd.LCBD = sd(LCBD),
            min.LCBD = min(LCBD),
            max.LCBD = max(LCBD))
sum2

#VegType  Avg.LCBD sd.LCBD  min.LCBD  max.LCBD
#1 Emergent   0.0494  0.0178   0.0189   0.0728
#2 Invaded    0.0459  0.0179   0.0319   0.0934
#3 Meadow     0.0548  0.0167   0.0319   0.0875

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


#### Temporal Beta Diversity Index ####

TBI1 <- TBI(spp2014, spp2015, method = "%difference", nperm = 999, test.t.perm = FALSE)
TBI1

TBI1$BCD.mat

compare2$TBI <- TBI1$BCD.mat
compare2$p.TBI <- TBI1$p.TBI
compare2

write.csv(compare2, "Data/LCBD_TBI_data_rel.csv")

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




