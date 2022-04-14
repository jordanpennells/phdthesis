

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



morfi_process_benchmarking <- read_csv("morfi_process_benchmarking.csv")


#PCA by factor all biomass

morfi_benchmarking_PCA <- subset(morfi_process_benchmarking,select = -c(Process,Level,Morfi.Rep))

PCA.benchmarking <- prcomp(~ ., data = morfi_benchmarking_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.benchmarking)

PCA.benchmarking

fviz_eig(PCA.benchmarking)

ggbiplot(PCA.benchmarking)


#Naming factors

PCA.benchmarking.process <- c(rep("HPH",12),rep("Bead Mill",8),rep("Silverson",20),rep("Extrusion-10%",12),rep("Extrusion-15%",12),rep("Extrusion-20%",12),rep("Extrusion-25%",12))
PCA.benchmarking.process <- factor(PCA.benchmarking.process, levels = c("HPH","Bead Mill","Silverson","Extrusion-10%","Extrusion-15%","Extrusion-20%","Extrusion-25%"))

PCA.benchmarking.level <- c(rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("H",4),rep("LL",4),rep("L",4),rep("M",4),rep("H",4),rep("HH",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4))
PCA.benchmarking.level <- factor(PCA.benchmarking.level, levels = c("LL","L","M","H","HH"))

pca_process <- ggbiplot(PCA.benchmarking,varname.size = 4,obs.scale = 0.4, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.process) + scale_color_manual(name = 'Process', values=c("darkgoldenrod1","darkgrey","chartreuse2","brown", "blue", "darkorange1", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_process + xlim(-3,6.2) + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


pca_HPH <- ggbiplot(PCA.benchmarking,varname.size = 4,obs.scale = 0.9, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.4, groups = PCA.benchmarking.level) + 
  scale_color_manual(name = 'Level', values=c("lightgreen","green", "orange", "red","darkred")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
  # geom_text(label = morfi_process_benchmarking$Level, alpha = 0.6)

pca_HPH + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))




#Hierarchical Clustering


morfi_process_benchmarking$biomass_ID <- with(morfi_process_benchmarking,paste(Process,Level, sep= "-"))

morfi_process_benchmarking_means <- morfi_process_benchmarking %>% group_by(biomass_ID) %>% summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
morfi_process_benchmarking_means <- morfi_process_benchmarking_means %>% dplyr:: select(1,3:23)

morfi_process_benchmarking_means <- as.data.frame(morfi_process_benchmarking_means)

rownames(morfi_process_benchmarking_means) <- morfi_process_benchmarking_means$biomass_ID
morfi_process_benchmarking_means <- morfi_process_benchmarking_means[,-1]

morfi_process_benchmarking_means <- scale(morfi_process_benchmarking_means,center=TRUE,scale=TRUE)

biomass.dist <- dist(morfi_process_benchmarking_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"ward.D2")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)



morfi.dist <- dist(t(morfi_process_benchmarking_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"ward.D2")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)


















