
# Libraries ---------------------------------------------------------------

library(vegan)
library(tidyverse)

library(mvabund)

library(car)
library(agricolae)

library(performance)
library(see)

library(lmPerm)

library(patchwork)

# Data --------------------------------------------------------------------

species <- read.csv("Data/Species matrix_raw.csv")

species$Year <- as.factor(species$Year)

species <- species %>% 
  mutate(Veg = fct_recode(VegType,
                              "Remnant" = "Emergent",
                              "Remnant" = "Meadow")) %>% 
  unite("VegYr", Year:VegType, remove = FALSE) %>% 
  unite("veg.year", Veg,Year, remove = FALSE) %>% 
  mutate(S = specnumber(species[,6:28]), # species richness
         Ab = rowSums(species[,6:28]), # total abundance
         H = diversity(species[,6:28]), # Shannon Index
         J = (H/log(specnumber(species[,6:28]))), # Pielou's evenness
         D = diversity(species[,6:28], "simpson"), #Simpson 1 - D
         I = diversity(species[,6:28], "invsimpson"), # Simpsons 1/D
         logS = log(S),
         logAb = log(Ab))
           
  
univariate <- species %>% 
  select(Site:Depth,Veg:logAb)

write.csv(univariate, "Data/species_univariate.csv")
write.csv(species, "Data/matrix_univariate.csv")


# ANOVAs ------------------------------------------------------------------

un.data <- read.csv("Data/species_univariate.csv")
un.data$Year <- as.factor(un.data$Year)

colnames(un.data)

## Histograms


par(mfrow = c(1,2))

hist(un.data$logAb, 
     xlab = "Total abundance", main = " ",
     border = "black",
     col = "white")

hist(un.data$logS, 
     xlab = "Total SpeciesRichness", main = " ",
     border = "black",
     col = "white")


## Abundance 

ab.mod <- lm(logAb ~ Veg * Year, data = un.data)
Anova(ab.mod, type = "3")

#Anova Table (Type III tests)

#Response: logAb
#            Sum Sq Df  F value    Pr(>F)    
#(Intercept) 45.357  1 444.6357 < 2.2e-16 ***
#Veg          0.250  1   2.4508  0.126218    
#Year         0.012  1   0.1162  0.735165    
#Veg:Year     0.967  1   9.4832  0.003956 ** 
#Residuals    3.672 36 

plot(ab.mod)

check_model(ab.mod) # not bad

## Species richness

s.mod <- lm(logS ~ Veg * Year, data = un.data)
Anova(s.mod, type = "3")


#Response: logS
#Sum Sq Df  F value  Pr(>F)    
#(Intercept) 17.6150  1 252.2236 < 2e-16 ***
#Veg          0.1462  1   2.0939 0.15654    
#Year         0.0228  1   0.3267 0.57117    
#Veg:Year     0.2739  1   3.9220 0.05534 .  
#Residuals    2.5142 36    

check_model(s.mod)

# Simpsons Inverse

i.mod <- lm(I ~ Veg * Year, data = un.data)
Anova(i.mod, type = "3")

#Response: I
#            Sum Sq Df  F value    Pr(>F)    
#(Intercept) 89.905  1 128.0930 2.057e-13 ***
#Veg          1.680  1   2.3940    0.1305    
#Year         0.343  1   0.4891    0.4888    
#Veg:Year     0.497  1   0.7083    0.4056    
#Residuals   25.267 36  

check_model(i.mod)

# Shannon-Weiner

h.mod <- lm(H ~ Veg * Year, data = un.data)
Anova(h.mod, type = "3")

#Response: H
#             Sum Sq Df  F value Pr(>F)    
#(Intercept) 13.4489  1 244.4013 <2e-16 ***
#Veg          0.1401  1   2.5461 0.1193    
#Year         0.0177  1   0.3216 0.5742    
#Veg:Year     0.0696  1   1.2645 0.2682    
#Residuals    1.9810 36 

check_model(h.mod)

## Pielous 

j.mod <- lm(J ~ Veg * Year, data = un.data)
Anova(j.mod, type = "3")

#Response: J
#             Sum Sq Df  F value Pr(>F)    
#(Intercept) 4.6166  1 745.9142 <2e-16 ***
#Veg         0.0005  1   0.0885 0.7677    
#Year        0.0000  1   0.0020 0.9644    
#Veg:Year    0.0141  1   2.2765 0.1401    
#Residuals   0.2228 36

check_model(j.mod)


# multi GLM ---------------------------------------------------------------

species <- read.csv("Data/matrix_univariate.csv")

spp <- species %>% 
  select(AMBI:AMRO)

env <- species %>% 
  select(Site:Depth,Veg)


spp.mv <- mvabund(spp)


spp.mod <- manyglm(spp.mv ~ Veg * factor(Year),
                   data = env, family = "poisson")

plot(spp.mod)


anova(spp.mod, p.uni = "adjusted") 

#Model: spp.mv ~ Veg * factor(Year)

#Multivariate test:
#                  Res.Df Df.diff   Dev Pr(>Dev)   
#(Intercept)          39                          
#Veg                  38       1 30.76    0.421   
#factor(Year)         37       1 70.01    0.003 **
#Veg:factor(Year)     36       1 33.92    0.028 * 


# Univariate Tests:
#                      AMBI           LEBI           MAWR              SORA            SWSP            VIRA            AMWO             COYE                RWBL
#                    Dev Pr(>Dev)    Dev Pr(>Dev)   Dev Pr(>Dev)      Dev Pr(>Dev)     Dev Pr(>Dev)    Dev Pr(>Dev)    Dev Pr(>Dev)     Dev Pr(>Dev)        Dev
#(Intercept)                                                                                                                                    
#Veg              2.853    0.879   1.427    0.988    0.555    0.997   0.016    0.997   0.215    0.997   0.946    0.994   2.408    0.970    1.117    0.994  6.308
#factor(Year)     1.046    0.997       0    1.000    2.289    0.962    0.34    0.997   0.038    1.000   1.328    0.997   1.386    0.997     0.91    0.997 33.931
#Veg:factor(Year)     0    0.980       0    0.996    0.736    0.964   3.819    0.603    0.49    0.964   2.969    0.678       0    0.964    0.102    0.964  9.162

#                                    YWAR           SOSP              CHSP             NOCA           WIFL            TRES           BARS           BANS
#                 Pr(>Dev)      Dev Pr(>Dev)      Dev Pr(>Dev)      Dev Pr(>Dev)     Dev Pr(>Dev)   Dev Pr(>Dev)    Dev Pr(>Dev)   Dev Pr(>Dev)   Dev
#(Intercept)                                                                                                                                      
#Veg                 0.450    0.047    0.997    0.062    0.997    0.713    0.997   0.713    0.997 0.713    0.994  0.447    0.997 3.841    0.742 1.427
#factor(Year)        0.002    1.007    0.997    1.359    0.997    1.386    0.997   1.386    0.997 1.386    0.997   0.04    1.000     0    1.000 2.772
#Veg:factor(Year)    0.141    0.486    0.964    0.734    0.964        0    0.980       0    0.980     0    0.980 10.918    0.098 1.726    0.852     0

#                                  EAKI           PUMA           CSWA           COGR          CLIFF           AMRO         
#                   Pr(>Dev)    Dev Pr(>Dev)   Dev Pr(>Dev)   Dev Pr(>Dev)   Dev Pr(>Dev)   Dev Pr(>Dev)   Dev Pr(>Dev)
#(Intercept)                                                                                                         
#Veg                 0.988  2.326    0.970 1.427    0.991 1.427    0.987 0.349    0.997 0.713    0.997 0.713    0.997
#factor(Year)        0.922 13.863    0.037 2.772    0.922     0    1.000     0    1.000 1.386    0.997 1.386    0.997
#Veg:factor(Year)    0.964      0    0.980     0    0.964     0    0.996 2.773    0.703     0    0.980     0    0.975

#Arguments:
#Test statistics calculated assuming uncorrelated response (for faster computation) 
#P-value calculated using 999 iterations via PIT-trap resampling.

# Figures -----------------------------------------------------------------

colour = c("Invaded" = "#fc8d62",
           "Remnant" = "#35978f")

shape = c("2014" = 21,
          "2015" = 24)

# Univariate figures ------------------------------------------------------

abundance <- ggplot(un.data,
                   aes(x = Veg, y = Ab)) + 
  geom_jitter(
    aes(shape = Year, fill =  Veg), 
    position = position_jitterdodge(jitter.width = 0.2,
                                    dodge.width = 0.6),
    size = 4.5,
    stroke = 1.5,
    alpha = 0.8) +
  theme_bw() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6),
    fill = "black") +
  labs(x = " ",
       y = expression(paste("Total Abundance"))) + 
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(0, 40)

richness <- ggplot(un.data,
                    aes(x = Veg, y = S)) + 
  geom_jitter(
    aes(shape = Year, fill =  Veg), 
    position = position_jitterdodge(jitter.width = 0.2,
                                    dodge.width = 0.6),
    size = 4.5,
    stroke = 1.5) +
  theme_classic() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6),
    fill = "black") +
  labs(x = " ",
       y = expression(paste("Species Richness"))) + 
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(0, 12)

simp <- ggplot(un.data,
               aes(x = Veg, y = I)) + 
  geom_jitter(
    aes(shape = Year, fill =  Veg), 
    position = position_jitterdodge(jitter.width = 0.2,
                                    dodge.width = 0.6),
    size = 4.5,
    stroke = 1.5,
    alpha = 0.8) +
  theme_bw() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6),
    fill = "black") +
  labs(x = " ",
       y = expression(paste("Simpson's Index (I/D)"))) + 
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(0,7)

shannon <- ggplot(un.data,
                    aes(x = Veg, y = H)) + 
  geom_jitter(
    aes(shape = Year, fill =  Veg), 
    position = position_jitterdodge(jitter.width = 0.2,
                                    dodge.width = 0.6),
    size = 4.5,
    stroke = 1.5,
    alpha = 0.8) +
  theme_bw() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6),
    fill = "black") +
  labs(x = " ",
       y = expression(paste("Shannon-Weiner Index"))) + 
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  theme(legend.position = "blank") +
  ylim(0,3)


pielou <- ggplot(un.data,
                 aes(x = Veg, y = J)) + 
  geom_jitter(
    aes(shape = Year, fill =  Veg), 
    position = position_jitterdodge(jitter.width = 0.2,
                                    dodge.width = 0.6),
    size = 4.5,
    stroke = 1.5,
    alpha = 0.8) +
  theme_bw() +
  stat_summary(
    aes(shape = Year),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6),
    fill = "black",
    show.legend = F) +
  labs(x = " ",
       y = expression(paste("Pielou's Evenness"))) + 
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "right") +
  guides(fill = "none") +
  ylim(0, 1)


panel <- abundance + richness + shannon + pielou +
  plot_annotation(tag_levels = 'A')

panel

ggsave("Figures/alpha_diversity_abundance.TIFF", panel)

figure2 <- year14 + year15 + 
  abundance + richness + shannon + pielou +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(nrow = 3)

ggsave("Figures/Figure2_iNEXT_uni.TIFF",
       figure3)


# Multivariate figures ----------------------------------------------------

spp.wide <- species %>% 
  relocate(Veg, .after = VegType) %>% 
  select(Site:AMRO)

spp.long <- spp.wide  %>% 
  pivot_longer(AMBI:AMRO, names_to = "Species", values_to = "count")

spp.long$Year <- as.factor(spp.long$Year)



cattail <- spp.long %>% 
  filter(VegType == "Emergent") %>% 
  filter(count > 0)


ctsum <- cattail %>% group_by(Year, Species) %>% 
  summarize(sum = sum(count))

ctsum$veg <- c("Emergent")
               
# total 
#1 AMBI        4
#2 BANS        1
#3 BARS        1
#4 COYE       49
#5 CSWA        1
#6 EAKI        4
#7 LEBI        2
#8 MAWR       51
#9 RWBL      133
#10 SWSP       50
#11 TRES        7
#12 VIRA        4
#13 WIFL        1
#14 YWAR        8


meadow <- spp.long %>% 
  filter(VegType == "Meadow")%>% 
  filter(count > 0)


meadow %>% group_by(Species) %>% 
  summarize(sum = sum(count))

#Species   sum
#<chr>   <int>
#1 AMRO        1
#2 BANS        1
#3 BARS        3
#4 CHSP        1
#5 CLIFF       1
#6 COGR        1
#7 COYE       33
#8 CSWA        1
#9 EAKI        5
#10 MAWR        6
#11 NOCA        1
#12 PUMA        2
#13 RWBL       84
#14 SORA        2
#15 SOSP        8
#16 SWSP       22
#17 TRES       12
#18 VIRA        2
#19 YWAR        9

mdsum <- meadow %>% group_by(Year, Species) %>% 
  summarize(sum = sum(count))

mdsum$veg <- c("Meadow")

phrag <- spp.long %>% 
  filter(VegType == "Invaded")%>% 
  filter(count > 0)


phrag %>% group_by(Species) %>% 
  summarize(sum = sum(count))

#Species   sum
#1 AMWO        1
#2 BARS        6
#3 COGR        1
#4 COYE       28
#5 EAKI        1
#6 MAWR       29
#7 RWBL       66
#8 SORA        1
#9 SOSP        4
#10 SWSP       34
#11 TRES        6
#12 VIRA        1
#13 YWAR        8

phr.sum <- phrag %>% group_by(Year, Species) %>% 
  summarize(sum = sum(count))

phr.sum$veg <- c("Phragmites")

test <- rbind(ctsum, mdsum)
birds.veg <- rbind(test, phr.sum)

birds.veg <- birds.veg %>% 
  pivot_wider(names_from = "Species",
              values_from = "sum") %>%
  replace(is.na(.), 0)

birds.veg <- birds.veg %>% t 

write.csv(birds.veg, "Data/birds_veg_appendix.csv")



ct <- ggplot(cattail, aes(x = Year, 
                         y = count,
                         shape = Year)) +
  facet_wrap(~Species) +
  ggtitle("Cattail") +
  geom_boxplot(lwd = 0.7,
               width = 0.5,
               position = position_dodge(0.8),
               show.legend = F) +
  geom_jitter(aes(stroke = 1,
                  shape = Year),
              size = 3,
              width = 0.2,
              fill = "#35978f") +
  labs(y = 'Count', x = ' ') +
  theme_bw() +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  guides(fill ="none")

md <- ggplot(meadow, aes(x = Year, 
                         y = count,
                         shape = Year)) +
  facet_wrap(~Species) +
  ggtitle("Meadow") +
  geom_boxplot(lwd = 0.7,
               width = 0.5,
               position = position_dodge(0.8),
               show.legend = F) +
  geom_jitter(aes(stroke = 1,
                  shape = Year),
              size = 3,
              width = 0.2,
              fill = "#35978f") +
  labs(y = 'Count', x = ' ') +
  theme_bw() +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  guides(fill ="none")
  

phr <- ggplot(phrag, aes(x = Year, y = count,
                         shape = Year)) +
  facet_wrap(~Species) +
  ggtitle("Phragmites") +
  geom_boxplot(lwd = 0.7,
               width = 0.5,
               position = position_dodge(0.8),
               show.legend = F) +
  geom_jitter(aes(stroke = 1,
                  fill = "#fc8d62"),
              size = 3,
              width = 0.2) +
  labs(y = 'Count', x = ' ') +
  theme_bw() +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill ="none")
  



veg <- ct + md + phr +
  plot_layout(ncol = 3)

veg


# Total count

ggplot(spp.long, aes(x = Veg, y = count + 1, 
                     fill = Veg,
                     shape = Year)) +
  geom_boxplot(lwd = 0.7,
               width = 0.5,
               position = position_dodge(0.8),
               show.legend = F) +
  geom_jitter(aes(stroke = 1),
              position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.8),
              size = 3) +
  labs(y = 'Abundance + 1', x = ' ') +
  theme_bw() +
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill ="none")



# Figure by species

sp.long <- spp.long %>% 
  filter(count > 0)



spp.facet <- ggplot(sp.long, aes(x = Veg, y = count, 
                                  fill = Veg,
                                  shape = Year)) +
  geom_boxplot(lwd = 0.7,
               width = 0.5,
               position = position_dodge(0.8),
               show.legend = F) +
  geom_jitter(aes(stroke = 1),
              position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.8),
              size = 4) +
  facet_wrap(~Species) +
  labs(y = 'Abundance', x = ' ') +
  theme_bw() +
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.75,0.07)) +
  guides(fill ="none")

spp.facet

ggsave("Figures/species_facet_wrap.TIFF",
       spp.facet)

