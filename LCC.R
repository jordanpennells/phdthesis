
library(factoextra)
library(ggplot2)
library(ggfortify)
library(devtools)
library(ggbiplot)
library(ggfortify)
library(ggrepel)
library(dplyr)
library(FactoMineR)
library(data.table)
library(stringr)
library(emmeans)
library(multcomp)
library(readr)
library(forcats)
library(tidyr)
library(gridExtra)
library(gclus)
library(ggdendro)
library(remotes)
library(tidyverse)
library(bibliometrix)
library(bibtex)
library(bib2df)
library(readxl)
library(janitor)
library(here)
library(quanteda)
library(collostructions)

#################################################################################### LCC


LCC <- read_csv("LCC_RStudio.csv")


#Energy v TI
ggplot(LCC) +
  aes(x = Energy, y = TI, shape = Study, color = Study, size = 10, fill = Study) +
  scale_shape_manual(values=c(15,16,17,18,25,19,11,24,8)) +
  geom_point() +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 14, face="bold"), legend.position = "right", legend.title=element_text(size=16,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), size = FALSE) +
  xlab("Energy (kWh/kg)") +
  ylab("Tensile Index (Nm/g)")



#Energy v TI - labels
ggplot(LCC) +
  aes(x = Energy, y = TI, shape = Study, color = Study, size = 10, fill = Study) +
  scale_shape_manual(values=c(15,16,17,18,25,19,11,24,8)) +
  geom_jitter(alpha = 0.2) +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 14, face="bold"), legend.position = "right", legend.title=element_text(size=16,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), size = FALSE) +
  xlab("Energy (kWh/kg)") +
  ylab("Tensile Index (Nm/g)") +
  geom_text_repel(label = LCC$Sample)



#Energy v TI - Ellipse
nomissing <- na.omit(LCC)
df <- nomissing
find_hull <- function(df) df[chull(df$TI, df$Energy), ]
hulls <- ddply(df, "Study", find_hull)


ggplot(LCC) +
  aes(x = Energy, y = TI, shape = Study, color = Study, size = 5, fill = Study) +
  scale_shape_manual(values=c(15,16,17,18,25,19,11,24,8)) +
  geom_point() +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 14, face="bold"), legend.position = "right", legend.title=element_text(size=16,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  guides(shape = guide_legend(override.aes = list(size = 1.5)), size = FALSE) +
  xlab("Energy (kWh/kg)") +
  ylab("Tensile Index (Nm/g)") +
  geom_polygon(data = hulls, alpha = 0.2, size = 1)







