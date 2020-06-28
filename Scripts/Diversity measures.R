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



####### Beta Diversity paritioning and TBI ######

library(adespatial)

Year2014 <- species %>% filter(Year == "2014")
Year2015 <- species %>% filter(Year == "2015")

spp2014 <- Year2014[ ,5:27]
env2014 <- Year2014[ ,1:4]

spp2015 <- Year2015[ ,5:27]
env2015 <- Year2015[ ,1:4]


##### Create Bray-Curtis dissimilarity matrix ###

spp2014bc <- vegdist(spp2014, method = "bray",
                     binary = FALSE)
spp2014bc


spp2015bc <- vegdist(spp2015, method = "bray",
                     binary = FALSE)
spp2015bc

## beta diversity 

BD2014 <- beta.div(spp2014bc, method = "percentdiff",
         sqrt.D = FALSE, samp = FALSE,
         nperm = 999)

Year2014$LCBD <- BD2014$LCBD

BD2014$LCBD[BD2014$LCBD > mean(BD2014$LCBD)]

#4           5           6           7           9          10         13         19 
#0.12124900  0.06212209  0.05502379  0.07485454  0.05221660 0.05759652 0.05109404 0.05577960 

YearF <- ggplot(Year2014, aes(y = LCBD, x = VegType)) +
  geom_point() +
  theme_classic()

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









#### TBI ####

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

# test difference that there is no sig. difference between times

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



