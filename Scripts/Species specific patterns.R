

# Packages ----------------------------------------------------------------

library(vegan)
library(tidyverse)
library(patchwork)

# Load data ---------------------------------------------------------------

species <- read.csv("Data/Species matrix_raw.csv")

species <- species %>% 
  mutate(VegType = fct_recode(VegType,
                              "Remnant" = "Emergent",
                              "Remnant" = "Meadow")) %>% 
  unite("VegYr", Year:VegType, remove = FALSE)


spp.long <- species %>% 
  pivot_longer(AMBI:AMRO, names_to = "Species", values_to = "count")


spp.long$Year <- as.factor(spp.long$Year)

colnames(spp.long)


colour = c("Invaded" = "#fc8d62",
           "Remnant" = "#35978f")

spp.facet <- ggplot(spp.long, aes(x = VegType, y = count + 1, 
                         fill = VegType,
                         shape = Year)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.6),
              size = 3) +
  facet_wrap(~Species) +
  scale_y_log10() +
  labs(y = 'Abundance + 1', x = ' ') +
  theme_bw() +
  scale_fill_manual(values = colour) +
  scale_shape_manual(values = c(21,24)) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.75,0.07)) +
  guides(fill ="none")

ggsave("Figures/species_facet_wrap.TIFF",
       spp.facet)



