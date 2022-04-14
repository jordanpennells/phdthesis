


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
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(csv)
library(corrplot)
library(RColorBrewer)

############################# MorFi Data Analysis ########################

morfi_benchmarking <- read_csv("morfi_benchmarking.csv")


#PCA by factor all biomass

morfi_benchmarking_PCA <- subset(morfi_benchmarking,select = -c(Biomass,Section,Level,Morfi.Rep))

PCA.benchmarking <- prcomp(~ ., data = morfi_benchmarking_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.benchmarking)

PCA.benchmarking

fviz_eig(PCA.benchmarking)

ggbiplot(PCA.benchmarking)


PCA.benchmarking.biomass <- c(rep("Banana",24),rep("Softwood",12),rep("Sugarcane",24),rep("Spinifex",12),rep("Sugargraze",48),rep("Yemen",48),rep("GreenleafBMR",36),rep("Graingrass",35))
PCA.benchmarking.biomass <- factor(PCA.benchmarking.biomass, levels = c("Banana","Softwood","Sugarcane","Spinifex","Sugargraze","Yemen","GreenleafBMR","Graingrass"))

PCA.benchmarking.HPH <- c(rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",3),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4))
PCA.benchmarking.HPH <- factor(PCA.benchmarking.HPH, levels = c("L","M","H"))


pca_biomass <- ggbiplot(PCA.benchmarking,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.biomass) + scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","darkgrey","chartreuse2","brown", "blue", "darkorange1", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


pca_HPH <- ggbiplot(PCA.benchmarking,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.HPH) + scale_color_manual(name = 'HPH Level', values=c("lightgreen", "orange", "red")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_HPH + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


#PCA by factor non-wood

morfi_benchmarking_nonwood <- read_csv("morfi_benchmarking_nonwood.csv")

morfi_benchmarking_nonwood_PCA <- subset(morfi_benchmarking_nonwood,select = -c(Biomass,Section,Level,Morfi.Rep))

PCA.benchmarking.nonwood <- prcomp(~ ., data = morfi_benchmarking_nonwood_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.benchmarking.nonwood)

PCA.benchmarking.nonwood

fviz_eig(PCA.benchmarking.nonwood)

ggbiplot(PCA.benchmarking.nonwood)


PCA.benchmarking.nonwood.biomass <- c(rep("Banana",24),rep("Sugarcane",24),rep("Spinifex",12),rep("Sugargraze",48),rep("Yemen",48),rep("GreenleafBMR",36),rep("Graingrass",35))
PCA.benchmarking.nonwood.biomass <- factor(PCA.benchmarking.nonwood.biomass, levels = c("Banana","Sugarcane","Spinifex","Sugargraze","Yemen","GreenleafBMR","Graingrass"))

PCA.benchmarking.nonwood.HPH <- c(rep("L",4),rep("M",4),rep("H",4),  #Banana-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #Banana-Stem
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugarcane-Bagasse
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugarcane-Mulch
                                  rep("L",4),rep("M",4),rep("H",4),  #Spinifex
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugargraze-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugargraze-Sheath
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugargraze-Stem<1m
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugargraze-Stem<1m
                                  rep("L",4),rep("M",4),rep("H",4),  #Yemen-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #Yemen-Sheath
                                  rep("L",4),rep("M",4),rep("H",4),  #Yemen-Stem<1m
                                  rep("L",4),rep("M",4),rep("H",4),  #Yemen-Stem<1m
                                  rep("L",4),rep("M",4),rep("H",4),  #GreenleafBMR-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #GreenleafBMR-Sheath
                                  rep("L",4),rep("M",4),rep("H",4),  #GreenleafBMR-Stem
                                  rep("L",4),rep("M",4),rep("H",4),  #GG-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #GG-Sheath
                                  rep("L",4),rep("M",4),rep("H",4))  #GG-Stem
PCA.benchmarking.nonwood.HPH <- factor(PCA.benchmarking.nonwood.HPH, levels = c("L","M","H"))


pca_biomass_nonwood <- ggbiplot(PCA.benchmarking.nonwood,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.biomass) + scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","chartreuse2","brown", "blue", "darkorange1", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nonwood + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


pca_HPH_nonwood <- ggbiplot(PCA.benchmarking.nonwood,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.HPH) + scale_color_manual(name = 'HPH Level', values=c("lightgreen", "orange", "red")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_HPH_nonwood + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


#PCA by factor non-wood + different sections

morfi_benchmarking_nonwood_sections <- read_csv("morfi_benchmarking_nonwood_sections.csv")

morfi_benchmarking_nonwood_sections_PCA <- subset(morfi_benchmarking_nonwood_sections,select = -c(Biomass,Section,Level,Morfi.Rep))

PCA.benchmarking.nonwood.sections <- prcomp(~ ., data = morfi_benchmarking_nonwood_sections_PCA, na.action=na.omit, center = TRUE, scale=TRUE)



PCA.benchmarking.nonwood.biomass.sections <- c(rep("Banana-Leaf",12),
                                               rep("Banana-Stem",12),
                                               rep("Sugarcane-Bagasse",12),
                                               rep("Sugarcane-Mulch",12),
                                               rep("Spinifex-Leaf",12),
                                               rep("Sorghum-Tall-Leaf",24),
                                               rep("Sorghum-Tall-Sheath",24),
                                               rep("Sorghum-Tall-Stem(<1m)",24),
                                               rep("Sorghum-Tall-Stem(>1m)",24),
                                               rep("Sorghum-Short-Leaf",23),
                                               rep("Sorghum-Short-Sheath",24),
                                               rep("Sorghum-Short-Stem",24))

PCA.benchmarking.nonwood.biomass.sections <- factor(PCA.benchmarking.nonwood.biomass.sections, levels = c("Banana-Leaf",
                                                                                                 "Banana-Stem",
                                                                                                 "Sugarcane-Bagasse",
                                                                                                 "Sugarcane-Mulch",
                                                                                                 "Spinifex-Leaf",
                                                                                                 "Sorghum-Tall-Leaf",
                                                                                                 "Sorghum-Tall-Sheath",
                                                                                                 "Sorghum-Tall-Stem(<1m)",
                                                                                                 "Sorghum-Tall-Stem(>1m)",
                                                                                                 "Sorghum-Short-Leaf",
                                                                                                 "Sorghum-Short-Sheath",
                                                                                                 "Sorghum-Short-Stem"))

PCA.benchmarking.nonwood.HPH.sections <- c(rep("L",4),rep("M",4),rep("H",4),  #Banana-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),  #Banana-Stem
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugarcane-Bagasse
                                  rep("L",4),rep("M",4),rep("H",4),  #Sugarcane-Mulch
                                  rep("L",4),rep("M",4),rep("H",4),  #Spinifex
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Tall-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Tall-Sheath
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Tall-Stem<1m
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Tall-Stem>1m
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Short-Leaf
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),  #Sorghum-Short-Sheath
                                  rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4))  #Sorghum-Short-Stem
                                

PCA.benchmarking.nonwood.HPH.sections <- factor(PCA.benchmarking.nonwood.HPH.sections, levels = c("L","M","H"))


pca_biomass_nonwood_sections <- ggbiplot(PCA.benchmarking.nonwood.sections,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.biomass.sections) + 
  scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","darkgoldenrod3","chartreuse3","chartreuse3","brown", "deepskyblue","dodgerblue2","dodgerblue4","darkblue", "orchid1", "darkorchid1", "darkorchid4")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) + 
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") + 
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nonwood_sections + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))





#### Morfi PCA - Low only #############################################################

morfi_benchmarking_nonwood_sections_Low <- filter(morfi_benchmarking_nonwood_sections, Level =="L")
morfi_benchmarking_nonwood_sections_Low$ID <- with(morfi_benchmarking_nonwood_sections_Low,paste(Biomass,Section, sep= "-"))
PCA.benchmarking.nonwood.sections_Low_Group <- morfi_benchmarking_nonwood_sections_Low$ID
morfi_benchmarking_nonwood_sections_Low_PCA <- subset(morfi_benchmarking_nonwood_sections_Low,select = -c(Biomass,Section,Level,Morfi.Rep,ID))

PCA.benchmarking.nonwood.sections_Low_Group <- factor(PCA.benchmarking.nonwood.sections_Low_Group, levels = c("Banana-Leaf",
                                                                                                          "Banana-Stem",
                                                                                                          "Sugarcane-Bagasse",
                                                                                                          "Sugarcane-Mulch",
                                                                                                          "Spinifex-Leaf",
                                                                                                          "Sugargraze-Leaf",
                                                                                                          "Sugargraze-Sheath",
                                                                                                          "Sugargraze-Stem(<1m)",
                                                                                                          "Sugargraze-Stem(>1m)",
                                                                                                          "Yemen-Leaf",
                                                                                                          "Yemen-Sheath",
                                                                                                          "Yemen-Stem(<1m)",
                                                                                                          "Yemen-Stem(>1m)",
                                                                                                          "GreenleafBMR-Leaf",
                                                                                                          "GreenleafBMR-Sheath",
                                                                                                          "GreenleafBMR-Stem",
                                                                                                          "Graingrass-Leaf",
                                                                                                          "Graingrass-Sheath",
                                                                                                          "Graingrass-Stem"))

PCA.benchmarking.nonwood.sections_Low <- prcomp(~ ., data = morfi_benchmarking_nonwood_sections_Low_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

pca_biomass_nonwood_sections <- ggbiplot(PCA.benchmarking.nonwood.sections_Low,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.sections_Low_Group) + 
  scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","darkgoldenrod3",
                                                "chartreuse3","chartreuse3",
                                                "brown", 
                                                "deepskyblue","deepskyblue1","deepskyblue2","deepskyblue3",
                                                "orangered","orangered2","orangered3","orangered4",
                                                "springgreen","springgreen2","springgreen4",
                                                "orchid1", "darkorchid1", "darkorchid4")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) + 
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") + 
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nonwood_sections + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



###################################################################################

#### Morfi PCA - Low only - SECTIONS!!! #############################################################

morfi_benchmarking_nonwood_sections_Low <- filter(morfi_benchmarking_nonwood_sections, Level =="L")
morfi_benchmarking_nonwood_sections_Low$ID <- with(morfi_benchmarking_nonwood_sections_Low,paste(Biomass,Section, sep= "-"))
PCA.benchmarking.nonwood.sections_Low_Group_Sections <- morfi_benchmarking_nonwood_sections_Low$Section
morfi_benchmarking_nonwood_sections_Low_PCA <- subset(morfi_benchmarking_nonwood_sections_Low,select = -c(Biomass,Section,Level,Morfi.Rep,ID))

PCA.benchmarking.nonwood.sections_Low_Group <- factor(PCA.benchmarking.nonwood.sections_Low_Group, levels = c("Banana-Leaf",
                                                                                                              "Banana-Stem",
                                                                                                              "Sugarcane-Bagasse",
                                                                                                              "Sugarcane-Mulch",
                                                                                                              "Spinifex-Leaf",
                                                                                                              "Sugargraze-Leaf",
                                                                                                              "Sugargraze-Sheath",
                                                                                                              "Sugargraze-Stem(<1m)",
                                                                                                              "Sugargraze-Stem(>1m)",
                                                                                                              "Yemen-Leaf",
                                                                                                              "Yemen-Sheath",
                                                                                                              "Yemen-Stem(<1m)",
                                                                                                              "Yemen-Stem(>1m)",
                                                                                                              "GreenleafBMR-Leaf",
                                                                                                              "GreenleafBMR-Sheath",
                                                                                                              "GreenleafBMR-Stem",
                                                                                                              "Graingrass-Leaf",
                                                                                                              "Graingrass-Sheath",
                                                                                                              "Graingrass-Stem"))

PCA.benchmarking.nonwood.sections_Low <- prcomp(~ ., data = morfi_benchmarking_nonwood_sections_Low_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

pca_biomass_nonwood_sections <- ggbiplot(PCA.benchmarking.nonwood.sections_Low,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.sections_Low_Group_Sections) + 
  scale_color_manual(name = 'Biomass', values=c("chartreuse3","springgreen4","chartreuse3",
                                                "brown", 
                                                "deepskyblue","deepskyblue1","deepskyblue2","deepskyblue3",
                                                "orangered","orangered2","orangered3","orangered4",
                                                "springgreen","springgreen2",
                                                "orchid1", "darkorchid1", "darkorchid4")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) + 
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") + 
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nonwood_sections + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



###################################################################################

#### Morfi PCA - High only #############################################################

morfi_benchmarking_nonwood_sections_High <- filter(morfi_benchmarking_nonwood_sections, Level =="H")
morfi_benchmarking_nonwood_sections_High$ID <- with(morfi_benchmarking_nonwood_sections_High,paste(Biomass,Section, sep= "-"))
PCA.benchmarking.nonwood.sections_High_Group <- morfi_benchmarking_nonwood_sections_High$ID
morfi_benchmarking_nonwood_sections_High_PCA <- subset(morfi_benchmarking_nonwood_sections_High,select = -c(Biomass,Section,Level,Morfi.Rep,ID))

PCA.benchmarking.nonwood.sections_High_Group <- factor(PCA.benchmarking.nonwood.sections_High_Group, levels = c("Banana-Leaf",
                                                                                                              "Banana-Stem",
                                                                                                              "Sugarcane-Bagasse",
                                                                                                              "Sugarcane-Mulch",
                                                                                                              "Spinifex-Leaf",
                                                                                                              "Sugargraze-Leaf",
                                                                                                              "Sugargraze-Sheath",
                                                                                                              "Sugargraze-Stem(<1m)",
                                                                                                              "Sugargraze-Stem(>1m)",
                                                                                                              "Yemen-Leaf",
                                                                                                              "Yemen-Sheath",
                                                                                                              "Yemen-Stem(<1m)",
                                                                                                              "Yemen-Stem(>1m)",
                                                                                                              "GreenleafBMR-Leaf",
                                                                                                              "GreenleafBMR-Sheath",
                                                                                                              "GreenleafBMR-Stem",
                                                                                                              "Graingrass-Leaf",
                                                                                                              "Graingrass-Sheath",
                                                                                                              "Graingrass-Stem"))

PCA.benchmarking.nonwood.sections_High <- prcomp(~ ., data = morfi_benchmarking_nonwood_sections_High_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

pca_biomass_nonwood_sections <- ggbiplot(PCA.benchmarking.nonwood.sections_High,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.sections_High_Group) + 
  scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","darkgoldenrod3",
                                                "chartreuse3","chartreuse3",
                                                "brown", 
                                                "deepskyblue","deepskyblue1","deepskyblue2","deepskyblue3",
                                                "orangered","orangered2","orangered3","orangered4",
                                                "springgreen","springgreen2","springgreen4",
                                                "orchid1", "darkorchid1", "darkorchid4")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) + 
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") + 
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nonwood_sections + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



###################################################################################







pca_HPH_nonwood_sections <- ggbiplot(PCA.benchmarking.nonwood.sections,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.nonwood.HPH.sections) + 
  scale_color_manual(name = 'HPH Level', values=c("lightgreen", "orange", "red")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) + 
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") + 
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_HPH_nonwood_sections + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


#Hierarchical Clustering - All Biomass

morfi_benchmarking$biomass_ID <- with(morfi_benchmarking,paste(Biomass,Section, sep= "-"))

morfi_benchmarking_means <- morfi_benchmarking %>% group_by(biomass_ID) %>% dplyr::summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
morfi_benchmarking_means <- morfi_benchmarking_means %>% dplyr:: select(1,3:23)

morfi_benchmarking_means <- as.data.frame(morfi_benchmarking_means)

rownames(morfi_benchmarking_means) <- morfi_benchmarking_means$biomass_ID
morfi_benchmarking_means <- morfi_benchmarking_means[,-1]

morfi_benchmarking_means <- scale(morfi_benchmarking_means,center=TRUE,scale=TRUE)



### Ward's Method

      #Biomass samples

biomass.dist <- dist(morfi_benchmarking_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"ward.D2")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)

      #MorFi parameters

morfi.dist <- dist(t(morfi_benchmarking_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"ward.D2")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)


### Single Linkage


biomass.dist <- dist(morfi_benchmarking_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"single")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)




morfi.dist <- dist(t(morfi_benchmarking_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"single")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)


### Average Linkage


biomass.dist <- dist(morfi_benchmarking_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"average")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)   + 
  labs(title = "Dendrogram - Biomass Benchmarking - Fibre Morphology")


############################### Average Linkage #################################

morfi.dist <- dist(t(morfi_benchmarking_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"average")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)


### Complete Linkage

biomass.dist <- dist(morfi_benchmarking_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"complete")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)




morfi.dist <- dist(t(morfi_benchmarking_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"complete")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)











############################# Nanopaper Data Analysis ########################


#PCA - Biomass - Nanopaper

nanopaper_benchmarking_biomass <- read_csv("nanopaper_benchmarking_biomass.csv")

nanopaper_benchmarking_biomass_PCA <- subset(nanopaper_benchmarking_biomass,select = -c(Biomass,Section,Level,Nanopaper,Strip))

PCA.benchmarking.biomass.nanopaper <- prcomp(~ ., data = nanopaper_benchmarking_biomass_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.benchmarking.biomass.nanopaper)

PCA.benchmarking.biomass.nanopaper

fviz_eig(PCA.benchmarking.biomass.nanopaper)

ggbiplot(PCA.benchmarking.biomass.nanopaper)


PCA.benchmarking.biomass.nanopaper.biomass <- c(rep("Banana-Leaf",43),
                                        rep("Banana-Stem",35),
                                        rep("Sugarcane-Bagasse",42),
                                        rep("Sugarcane-Mulch",43),
                                        rep("Spinifex-Leaf",48),
                                        rep("Sugargraze",166),
                                        rep("Yemen",179),
                                        rep("GreenleafBMR",110),
                                        rep("Graingrass",127))
PCA.benchmarking.biomass.nanopaper.biomass <- factor(PCA.benchmarking.biomass.nanopaper.biomass, levels = c("Banana-Leaf","Banana-Stem","Softwood","Sugarcane-Bagasse","Sugarcane-Mulch","Spinifex-Leaf","Sugargraze","Yemen","GreenleafBMR","Graingrass"))

PCA.benchmarking.HPH.nanopaper <- c(rep("L",15),rep("M",14),rep("H",14),
                                    rep("L",13),rep("M",15),rep("H",7),
                                    rep("L",16),rep("M",15),rep("H",11),
                                    rep("L",16),rep("M",13),rep("H",14),
                                    rep("L",16),rep("M",16),rep("H",16),
                                    rep("L",7),rep("M",8),rep("H",12),rep("L",16),rep("M",15),rep("H",16),rep("L",16),rep("M",16),rep("H",15),rep("L",15),rep("M",15),rep("H",15),
                                    rep("L",13),rep("M",16),rep("H",11),
                                    rep("L",16),rep("M",16),rep("H",15),
                                    rep("L",16),rep("M",16),rep("H",16),
                                    rep("L",15),rep("M",14),rep("H",15),
                                    rep("L",8),rep("M",12),rep("H",11),
                                    rep("L",13),rep("M",13),rep("H",12),
                                    rep("L",14),rep("M",16),rep("H",11),
                                    rep("L",14),rep("M",14),rep("H",14),
                                    rep("L",12),rep("M",15),rep("H",14),
                                    rep("L",15),rep("M",14),rep("H",15))
PCA.benchmarking.HPH.nanopaper <- factor(PCA.benchmarking.HPH.nanopaper, levels = c("L","M","H"))


pca_biomass_nanopaper <- ggbiplot(PCA.benchmarking.biomass.nanopaper,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.biomass.nanopaper.biomass) + 
  scale_color_manual(name = 'Biomass', values=c("darkgrey","black","chartreuse2","chartreuse4","brown", "blue", "darkorange1", "darkgreen","purple")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) +
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") +
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_biomass_nanopaper + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


pca_HPH_nanopaper <- ggbiplot(PCA.benchmarking.biomass.nanopaper,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.benchmarking.HPH.nanopaper) + scale_color_manual(name = 'HPH Level', values=c("lightgreen", "orange", "red")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
pca_HPH_nanopaper + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


#Nanopaper Hierarchical clustering





nanopaper_benchmarking_biomass_PCA_means <- as.data.frame(nanopaper_benchmarking_biomass_PCA_means)

rownames(nanopaper_benchmarking_biomass_PCA_means) <- nanopaper_benchmarking_biomass_PCA_means$biomass_ID
nanopaper_benchmarking_biomass_PCA_means <- nanopaper_benchmarking_biomass_PCA_means[,-1]

nanopaper_benchmarking_biomass_PCA_means <- scale(nanopaper_benchmarking_biomass_PCA_means,center=TRUE,scale=TRUE)

nanopaper.dist <- dist(t(nanopaper_benchmarking_biomass_PCA_means),"euclidean")

nanopaper.hc <- hclust(nanopaper.dist,"ward.D2")

nanopaper.hc.opt <- reorder.hclust(nanopaper.hc,nanopaper.dist)

ggdendrogram(nanopaper.hc.opt,rotate = TRUE)













############################# MorFi vs Nanopaper Data Analysis ########################

#Correlation Table - Biomass - Morfi vs Nanopaper

morfi_benchmarking_nonwood <- read_csv("morfi_benchmarking_nonwood.csv")

morfi_benchmarking_nonwood$biomass_ID <- with(morfi_benchmarking_nonwood,paste(Biomass,Section,Level, sep= "-"))
morfi_benchmarking_nonwood <- morfi_benchmarking_nonwood[,c(ncol(morfi_benchmarking_nonwood),1:(ncol(morfi_benchmarking_nonwood)-1))]

morfi_benchmarking_nonwood_PCA <- subset(morfi_benchmarking_nonwood,select = -c(Biomass,Section,Level,Morfi.Rep))

morfi_benchmarking_nonwood_PCA_means <- morfi_benchmarking_nonwood_PCA %>% group_by(biomass_ID) %>% dplyr::summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
morfi_benchmarking_nonwood_PCA_means <- morfi_benchmarking_nonwood_PCA_means %>% dplyr:: select(1:22)

morfi_benchmarking_nonwood_PCA_means <- as.data.frame(morfi_benchmarking_nonwood_PCA_means)

rownames(morfi_benchmarking_nonwood_PCA_means) <- morfi_benchmarking_nonwood_PCA_means$biomass_ID
morfi_benchmarking_nonwood_PCA_means <- morfi_benchmarking_nonwood_PCA_means[,-1]

morfi_benchmarking_nonwood_PCA_means <- scale(morfi_benchmarking_nonwood_PCA_means)

morfi_df <- morfi_benchmarking_nonwood_PCA_means


nanopaper_benchmarking_biomass <- read_csv("nanopaper_benchmarking_biomass.csv")

nanopaper_benchmarking_biomass$biomass_ID <- with(nanopaper_benchmarking_biomass,paste(Biomass,Section,Level, sep= "-"))
nanopaper_benchmarking_biomass <- nanopaper_benchmarking_biomass[,c(ncol(nanopaper_benchmarking_biomass),1:(ncol(nanopaper_benchmarking_biomass)-1))]

nanopaper_benchmarking_biomass_PCA <- subset(nanopaper_benchmarking_biomass,select = -c(Biomass,Section,Level,Nanopaper,Strip))

nanopaper_benchmarking_biomass_PCA_means <- nanopaper_benchmarking_biomass_PCA %>% group_by(biomass_ID) %>% dplyr::summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
nanopaper_benchmarking_biomass_PCA_means <- nanopaper_benchmarking_biomass_PCA_means %>% dplyr:: select(1:11)

nanopaper_benchmarking_biomass_PCA_means <- as.data.frame(nanopaper_benchmarking_biomass_PCA_means)

rownames(nanopaper_benchmarking_biomass_PCA_means) <- nanopaper_benchmarking_biomass_PCA_means$biomass_ID
nanopaper_benchmarking_biomass_PCA_means <- nanopaper_benchmarking_biomass_PCA_means[,-1]

nanopaper_benchmarking_biomass_PCA_means <- scale(nanopaper_benchmarking_biomass_PCA_means)

nanopaper_df <- nanopaper_benchmarking_biomass_PCA_means


cor(morfi_benchmarking_nonwood_PCA_means,nanopaper_benchmarking_biomass_PCA_means)


# corrr::correlate(bind_cols(morfi_benchmarking_nonwood_PCA_means[, -1], nanopaper_benchmarking_biomass_PCA_means[, -1])) %>%
#   filter(rowname %in% colnames(morfi_benchmarking_nonwood_PCA_means)) %>%
#   select(one_of(c("rowname", colnames(nanopaper_benchmarking_biomass_PCA_means)[-1]))) %>%
#   as.data.frame() %>%
#   column_to_rownames("rowname") %>%
#   as.matrix() %>%
#   corrplot::corrplot(is.corr = FALSE)

  
benchmarking_biomass_morfi_nanopaper_corr <- cor(morfi_benchmarking_nonwood_PCA_means[,1:21],nanopaper_benchmarking_biomass_PCA_means[,6:10])  
  
write.csv(benchmarking_biomass_morfi_nanopaper_corr,"benchmarking_biomass_morfi_nanopaper_corr.csv", row.names = TRUE)  
  

# Corrplot MorFi

morfi_benchmarking_nonwood_PCA_means <- as.data.frame(morfi_benchmarking_nonwood_PCA_means)

morfi_benchmarking_nonwood_PCA_means <- morfi_benchmarking_nonwood_PCA_means %>% dplyr::select(1:21)

colnames(morfi_benchmarking_nonwood_PCA_means)[8] <- "fibre_coarse"
 
M = cor(morfi_benchmarking_nonwood_PCA_means,use="complete.obs")  

# corrplot(M, order = 'FPC', type = 'lower', col=c("red","darkorange","darkgoldenrod1","darkolivegreen1","darkolivegreen2", "chartreuse4"))  
  
corrplot(M, order = 'FPC', 
         type = 'lower', 
         tl.cex = 1.2,
         col=brewer.pal(n=8, name="RdYlGn"), 
         tl.col="black", 
         diag = TRUE,
         mar=c(1,1,1,1))


corrplot(M, 
         method = 'color', 
         addCoef.col = 'white', 
         diag = TRUE, 
         order = 'hclust', 
         col=colorRampPalette(c("#a50026","#f46d43", "white","white","white","white","#66bd63","#006837"))(8), 
         addrect = 4, 
         tl.col="black",
         tl.cex = 1.2,
         rect.col = 'black', 
         rect.lwd = 5, 
         tl.pos = 'bt',
         number.cex=0.7,
         cl.pos = "b",
         cl.cex = 1.2,
         cl.ratio = 0.1,
         mar=c(0,0,1,0),
         insig='blank',
         addgrid.col = NA)






############################# MorFi-Nanopaper Double Dial Diagram ########################

   ### "https://uc-r.github.io/hc_clustering" test

df <- nanopaper_benchmarking_biomass_PCA_means
df <- na.omit(df)

df <- scale(df)

d <- dist(df, method = "euclidean")

hc2 <- agnes(df, method = "complete")
hc2$ac

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

ac <- function(x) {
  agnes(df, method = x)$ac
}

map_dbl(m, ac)

hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")


hc4 <- diana(nanopaper_benchmarking_biomass_PCA_means)
hc4$dc

pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")

hc5 <- hclust(d, method = "ward.D2")
sub_grp <- cutree(hc5, k = 4)
table(sub_grp)

nanopaper_benchmarking_biomass_PCA_means %>%
  mutate(cluster = sub_grp) %>%
  head

plot(hc5, cex = 0.8)
rect.hclust(hc5, k = 4, border = 2:5)

fviz_cluster(list(data = df, cluster = sub_grp))

hc_a <- agnes(df, method = "ward")
cutree(as.hclust(hc_a), k = 4)

hc_d <- diana(df)
cutree(as.hclust(hc_d), k = 4)



res.dist <- dist(df, method = "euclidean")

hc1 <- hclust(res.dist, method = "complete")
hc2 <- hclust(res.dist, method = "ward.D2")

dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2)


dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)

fviz_nbclust(df, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


### Real Diagram

morfi.dist <- dist(morfi_df,"euclidean")
nanopaper.dist <- dist(nanopaper_df, method = "euclidean")

hc1 <- hclust(morfi.dist, method = "ward.D2")
hc2 <- hclust(nanopaper.dist, method = "ward.D2")

plot(hc1, cex = 0.8)
plot(hc2, cex = 0.8)

dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)














############ Biomass Sustainability Triangle ##############

nanopaper_benchmarking_biomass <- read_csv("nanopaper_benchmarking_biomass.csv")


#Organising Tensile Index (TI) parameters


nanopaper_benchmarking_esr.TI <- dplyr::select(nanopaper_benchmarking_biomass,c(1:3,15))

nanopaper_benchmarking_esr.TI <- na.omit(nanopaper_benchmarking_esr.TI)

nanopaper_benchmarking_esr.TI$Level <- gsub("L",0,gsub("M",0.5,gsub("H",1,nanopaper_benchmarking_esr.TI$Level)))

nanopaper_benchmarking_esr.TI$ID = str_c(nanopaper_benchmarking_esr.TI$Biomass,'-',nanopaper_benchmarking_esr.TI$Section)

nanopaper_benchmarking_norm.TI <- nanopaper_benchmarking_esr.TI %>% group_by(ID,Level) %>% 
  dplyr::summarise(mean(TI_star))

nanopaper_benchmarking_esr.TI <- left_join(nanopaper_benchmarking_esr.TI,nanopaper_benchmarking_norm.TI,by=c("ID","Level"))

nanopaper_benchmarking_esr.TI$TI_norm <- with(nanopaper_benchmarking_esr.TI,(TI_star-min(`mean(TI_star)`))/(diff(range(`mean(TI_star)`))))


##ESR Plotting - Tensile Index

ggplot(nanopaper_benchmarking_esr.TI) +
  aes(x = Level, y = TI_norm, color = Biomass, shape = Section) +
  scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1","chartreuse2","brown", "blue", "darkorange1", "darkgreen","purple")) +
  scale_shape_manual(values=c(15,17,16,1,10,3,4)) + 
  theme_bw() + 
  geom_point() + 
  stat_summary(aes(group=ID), fun=mean, geom="line", colour="green") +
  ylab("Tensile Index (normalised)")



fits <- lmList(TI_norm ~ as.numeric(Level) | ID, data=nanopaper_benchmarking_esr.TI)

fits.coef <- coef(fits)

colnames(fits.coef) <- c("intercept","slope")

fits.coef$ID <- rownames(fits.coef)
IDsplit <- str_split(fits.coef$ID,"-",simplify = TRUE)
fits.coef$Biomass <- IDsplit[,1]
fits.coef$Section <- IDsplit[,2]

fits.coef$angle <- atan(fits.coef$intercept/fits.coef$slope)
fits.coef$mag <- sqrt((fits.coef$intercept^2)+(fits.coef$slope^2))

fits.coef[fits.coef<0] <- 0.005


fits.coef <- rbind(c(0.005,0.005,"Softwood-Stem","Softwood","Stem",0.005,0.005),fits.coef)

fits.coef$intercept <- as.numeric(fits.coef$intercept)
fits.coef$slope <- as.numeric(fits.coef$slope)

max(nanopaper_benchmarking_esr.TI$`mean(TI_star)`)

ggplot(fits.coef) +
  aes(x=slope,y=intercept, color = Biomass, shape = Section) + 
  scale_color_manual(name = 'Biomass', values=c("darkgoldenrod1",
                                                "purple",
                                                "darkgreen",
                                                "red",
                                                "brown", 
                                                "chartreuse2", 
                                                "blue",
                                                "darkorange1")) +
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1, size = 1) +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size=12), axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.4), plot.subtitle = element_text(hjust = 0.4, size = 14), plot.caption = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", color = "Variety", shape = "Section", title = "Processing Sustainability Triangle", subtitle = "Biomass Benchmarking - Tensile Index") +
  guides(color = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.5), shape = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.4))






################################################### TGA ##########################


Biomass_Benchmarking_TGA <-read_csv("Biomass_Benchmarking_TGA.csv")

Biomass_Benchmarking_TGA_Symbols <- pivot_longer(Biomass_Benchmarking_TGA,
                                        Step:Residue,
                                        names_to = "Metric",
                                        values_to = "Value")



p <- ggplot(subset(Biomass_Benchmarking_TGA_Symbols, Metric %in% c("Step"))) +
  aes(x=Biomass, y=Value) + 
  geom_point(size=5, aes(shape = Section, color = Metric)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  
  scale_color_manual(values=c("darkgreen", "orange2", "brown4","bisque3","grey50","black","black")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=-8, 
                   yend=0), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Thermogravimetric analysis", 
       subtitle="Biomass Benchmarking Study",
       xlab="Step (%)") +  
  coord_flip()

Step <- p + facet_grid( cols = vars(Metric), labeller = labeller(Metric = supp.labs)) +
  theme_bw() +
  guides(color = "none", shape = "none")






p <- ggplot(subset(Biomass_Benchmarking_TGA_Symbols, Metric %in% c("Onset"))) +
  aes(x=Biomass, y=Value) + 
  geom_point(size=5, aes(shape = Section, color = Metric)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("orange2", "brown4","bisque3","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=125, 
                   yend=max(Value)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title=" ", 
       subtitle=" ") +  
  coord_flip()

Onset <- p + facet_grid(cols = vars(Metric), labeller = labeller(Metric = supp.labs)) +
  theme_bw()  +
  theme(strip.text.y = element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(color = "none", shape = "none")





p <- ggplot(subset(Biomass_Benchmarking_TGA_Symbols, Metric %in% c("Residue"))) +
  aes(x=Biomass, y=Value) + 
  geom_point(size=5, aes(shape = Section, color = Metric)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("brown4","bisque3","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=22.5, 
                   yend=max(Value)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title=" ", 
       subtitle=" ") +  
  coord_flip()

supp.labs <- c("Step (%)", "Onset (°C)", "Residue (%)")
names(supp.labs) <- c("Step","Onset","Residue")

Residue <- p + facet_grid( cols = vars(Metric), labeller = labeller(Metric = supp.labs)) +
  theme_bw() +
  theme(strip.text.y = element_blank(), 
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
  guides(color = "none")




grid.arrange(Step, Onset, Residue, nrow=1)





################################################### Crystallinity ##########################


crystallinity <-read_csv("crystallinity.csv")

# crystallinity_benchmarking <- pivot_longer(crystallinity,
#                                                  Step:Residue,
#                                                  names_to = "Crystallinity",
#                                                  values_to = "pc")



# Biochem_Biomass_Symbols$Treatment = factor(Biochem_Biomass_Symbols$Treatment, 
#                                            levels=c("Raw","DL"), 
#                                            labels=c("Raw", "DL"))



p <- ggplot(subset(crystallinity)) +
  aes(x=Biomass, y=Crystallinity) + 
  geom_point(size=5, aes(shape = Type, color = Biomass)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("goldenrod1",
                              "darkorchid4",
                              "darkgreen",
                              "grey50",
                              "brown",
                              "chartreuse2",
                              "blue",
                              "orange2")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Crystallinity)), 
               linetype="dashed", 
               size=0.1) +
  ylab("Crystallinity(%)") +
  labs(title="Crystallinity Index", 
       subtitle="Biomass Benchmarking Study") +  
  coord_flip()

p + theme_bw() 


#   facet_grid(rows = vars(Treatment), cols = vars(Biochemical))














