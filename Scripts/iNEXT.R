
# library -----------------------------------------------------------------

library(iNEXT)
library(tidyverse)
library(patchwork)


# data --------------------------------------------------------------------
richness <- read.csv("Data/species_richness.csv")

remnant <- richness %>% 
  select(remnant2014, remnant2015)

invaded <- richness %>% 
  select(invaded2014, invaded2015)

v14 <- richness %>% 
  select(remnant2014, invaded2014)

v15 <- richness %>% 
  select(remnant2015, invaded2015)

v15 <- v15 %>% 
  rename(Remnant = remnant2015,
         Invaded = invaded2015)


# By year -----------------------------------------------------------------

y14 <- iNEXT(v14, q = c(0, 1,2), datatype = "abundance")
y15 <- iNEXT(v15, q = c(0, 1,2), datatype = "abundance")

year14 <- ggiNEXT(y14, type = 1, se = TRUE,
                  facet.var = "order",
                  color.var = "site") +
  ggtitle("2014") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") 

year15 <- ggiNEXT(y15, type = 1, se = TRUE, facet.var = "order")+
  ggtitle("2015") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") 
  
  
year14 + year15

# Remnant vs invaded ------------------------------------------------------

# remnant vegetation

t <- iNEXT(remnant, q = c(0, 1,2), datatype = "abundance")
ChaoRichness(remnant)

#             Observed Estimator Est_s.e. 95% Lower 95% Upper
#remnant2014       14    31.902   23.496    16.518   141.283
#remnant2015       20    28.142    8.253    21.567    62.301

ChaoSimpson(remnant)

#Observed Estimator Est_s.e. 95% Lower 95% Upper
#remnant2014    0.792     0.796    0.014     0.792     0.824
#remnant2015    0.727     0.729    0.022     0.727     0.773

ChaoShannon(remnant)

#          Observed Estimator Est_s.e 95% Lower 95% Upper
#remnant2014    1.828     1.904   0.096     1.828     2.092
#remnant2015    1.807     1.852   0.078     1.807     2.005

t.plot <- ggiNEXT(t, type = 1, se = TRUE, facet.var = "order")+
  ggtitle("Remnant Vegetation") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") 

# invaded vegetation

a <- iNEXT(invaded, q = c(0,1,2), datatype = "abundance")

ChaoRichness(invaded)

#            Observed Estimator Est_s.e. 95% Lower 95% Upper
#invaded2014       12    19.918   11.545    12.977    76.175
#invaded2015        8     8.247    0.723     8.013    12.694

ChaoSimpson(invaded)

#            Observed Estimator Est_s.e. 95% Lower 95% Upper
#invaded2014    0.809     0.817    0.022     0.809     0.861
#invaded2015    0.762     0.770    0.026     0.762     0.821

ChaoShannon(invaded)
#Observed Estimator Est_s.e 95% Lower 95% Upper
#invaded2014    1.927     2.026   0.127     1.927     2.274
#invaded2015    1.637     1.681   0.098     1.637     1.872

b.plot <- ggiNEXT(a, type = 1, se = TRUE, facet.var = "order") +
  ggtitle("Invaded vegetation") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") +
  ylim(0, 30) +
  xlim(0, 220)
  

plots <- b.plot + t.plot

ggsave("Figures/iNEXT_plots.jpeg",
       plots)



# by invaded, meadow, emergent --------------------------------------------

habrich <- read.csv("Data/species_richness_habitat.csv")
colnames(habrich)

meadow <- habrich %>% 
  select(marsh14, marsh15)

meadow  <- meadow  %>% 
  rename("2014" = "marsh14",
         "2015" = "marsh15")


emerg <- habrich %>% 
  select(emg14, emg15)


m <- iNEXT(meadow, q = c(0, 1,2), datatype = "abundance")

ChaoRichness(meadow)

#Observed Estimator Est_s.e. 95% Lower 95% Upper
#marsh14        8     8.982    1.843     8.089    18.905
#marsh15       17    21.136    4.851    17.668    42.604

ChaoSimpson(meadow)

#Observed Estimator Est_s.e. 95% Lower 95% Upper
#marsh14    0.749     0.763    0.035     0.749     0.831
#marsh15    0.759     0.765    0.034     0.759     0.832

ChaoShannon(meadow)

#        Observed Estimator Est_s.e 95% Lower 95% Upper
#marsh14    1.608     1.687   0.115     1.608     1.913
#marsh15    1.976     2.058   0.129     1.976     2.311


e <- iNEXT(emerg, q = c(0, 1,2), datatype = "abundance")

ChaoRichness(emerg)

#Observed Estimator Est_s.e. 95% Lower 95% Upper
#emg14       11    16.953    7.104    11.941    48.646
#emg15       12    19.958   11.602    12.982    76.488

ChaoShannon(emerg)

#Observed Estimator Est_s.e 95% Lower 95% Upper
#emg14    1.795     1.862   0.085     1.795     2.029
#emg15    1.542     1.593   0.091     1.542     1.773

ChaoSimpson(emerg)

#Observed Estimator Est_s.e. 95% Lower 95% Upper
#emg14    0.797     0.804    0.015     0.797     0.833
#emg15    0.693     0.697    0.026     0.693     0.747



m.plot <- ggiNEXT(m, type = 1, se = TRUE, facet.var = "order")+
  ggtitle("Meadow marsh") +
theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") +
  ylim(0, 30)


e.plot <- ggiNEXT(e, type = 1, se = TRUE, facet.var = "order")+
  ggtitle("Emergent marsh") +
theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#fc8d62","#35978f")) +
  scale_fill_manual(values = c("#fc8d62","#35978f")) +
  scale_shape_manual(values = c(18, 15)) +
  guides(fill ="none") +
  ylim(0, 30)

hab.plot <- e.plot + m.plot +
  plot_annotation(tag_levels = "A")

ggsave("Figures/habitat_iNEXT.jpeg")


hab.plot <- b.plot + t.plot + e.plot + m.plot + 
  plot_annotation(tag_levels = "A")

ggsave("Figures/habitat_iNEXT.TIFF",
       hab.plot,
       height = 10,
       width = 16)


citation("iNEXT")


# species accumulation curve ----------------------------------------------


# R SAC vignette ----------------------------------------------------------


??specaccum
# from vignette 

data(BCI)
sp1 <- specaccum(BCI) # exact method - Chiarucci et al 2008, mao tau estimate (Colwell et al. 2012)
sp2 <- specaccum(BCI, "random") # classic method (Gotelli and Colwell 2001)
sp2
summary(sp2)

plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
fitted(mod1)
plot(sp1)

## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)

## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)




# SAC w/ bir ddata --------------------------------------------------------

species <- read.csv("Data/Species matrix_raw.csv")
spp <- species %>% select(AMBI:AMRO)

sp1 <- specaccum(spp) # exact method - Chiarucci et al 2008, mao tau estimate (Colwell et al. 2012)
sp2 <- specaccum(spp, "random") # classic method (Gotelli and Colwell 2001)
sp2
summary(sp2)

plot(sp1, ci.type="poly", col="black", lwd=5, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="white", add=TRUE)

specslope(sp1, 20)

# at 20 sites
# 0.3277401

# at 30 sites
# 0.2274145
