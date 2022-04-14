

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
library(dendextend)
library(corrplot)
library(gridExtra)


################### Writing .csv data files ###################

sedimentation <- read_csv("sedimentation.csv")

sedimentation$ID <- str_c(sedimentation$Variety,'-',
                          sedimentation$Section,'-',
                          sedimentation$HPH)

sedimentation_PCA <- subset(sedimentation,select = -c(Sample,
                                                      Replicate,
                                                      Variety,
                                                      Section,
                                                      HPH))
WRV <- read_csv("WRV.csv")

WRV$ID <- str_c(WRV$Variety,'-',WRV$Section,'-',WRV$HPH)

WRV_PCA <- subset(WRV,select = -c(Sample,Variety,Section,HPH))

WRV_Sedi <- merge(x = WRV_PCA, y = sedimentation_PCA, by = "ID", all.x = TRUE)


nanopaper_CH8 <- read_csv("nanopaper_CH8.csv")

nanopaper_CH8$ID <- str_c(nanopaper_CH8$Variety,'-',nanopaper_CH8$Section,'-',nanopaper_CH8$HPH)







################### Cleaning data for PCA ###################


Nanopaper_WRV_Sedi <- merge(x = nanopaper_CH8, y = WRV_Sedi, by = "ID", all.x = TRUE)

Nanopaper_WRV_Sedi_na <- Nanopaper_WRV_Sedi[complete.cases(Nanopaper_WRV_Sedi), ] 

Nanopaper_WRV_Sedi2 <- subset(Nanopaper_WRV_Sedi,select = -c(ID,Sample,Variety,Section,HPH,nanopaper.rep))







################### PCA: Nanopaper + Fibre-Water ###################


Strength.WR_PCA <- prcomp(~ ., data=Nanopaper_WRV_Sedi2, na.action=na.omit, center = TRUE, scale=TRUE)


summary(Strength.WR_PCA)

fviz_eig(Strength.WR_PCA)

ggbiplot(Strength.WR_PCA)


Strength.WR.variety <- c(rep("Graingrass",144),rep("GreenleafBMR",136),rep("Sugargraze",176),rep("Yemen",192))
Strength.WR.variety <- factor(Strength.WR.variety, levels = c("Sugargraze","Yemen","GreenleafBMR","Graingrass"))

# p_variety <- ggbiplot(Strength.WR_PCA,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = Strength.WR.variety) + scale_color_manual(name = 'Varieties', values=c("blue", "orange", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
# p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



Strength.WR.section <- Nanopaper_WRV_Sedi_na$Section

p_section <- ggbiplot(Strength.WR_PCA,
                      varname.size = 3,
                      obs.scale = 0.5, 
                      ellipse=TRUE, 
                      ellipse.prob = 0.68, 
                      alpha = 0.5, 
                      groups = Strength.WR.section) + 
  scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_section + theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill="white",colour = "black", size = 2)) +
  ggtitle("PCA: Nanopaper + Fibre-Water Interaction Metrics")













################### PCA: MorFi + Fibre-Water ###################


morfi_cnf_CH8 <- read_csv("morfi_cnf.csv")

morfi_cnf_CH8$ID <- str_c(morfi_cnf_CH8$Variety,'-',morfi_cnf_CH8$Section,'-',morfi_cnf_CH8$HPH)

morfi_WRV_Sedi <- merge(x = morfi_cnf_CH8, y = WRV_Sedi, by = "ID", all.x = TRUE)

morfi_WRV_Sedi_na <- morfi_WRV_Sedi[complete.cases(morfi_WRV_Sedi), ] 

morfi_WRV_Sedi2 <- subset(morfi_WRV_Sedi,select = -c(ID,Sample_name,Variety,Replicate,Section,HPH,Morfi.Rep))



morfi_WRV_Sedi2_PCA <- prcomp(~ ., data=morfi_WRV_Sedi2, na.action=na.omit, center = TRUE, scale=TRUE)



summary(morfi_WRV_Sedi2_PCA)

fviz_eig(morfi_WRV_Sedi2_PCA)

ggbiplot(morfi_WRV_Sedi2_PCA)


# Strength.WR.variety <- c(rep("Graingrass",144),rep("GreenleafBMR",136),rep("Sugargraze",176),rep("Yemen",192))
# Strength.WR.variety <- factor(Strength.WR.variety, levels = c("Sugargraze","Yemen","GreenleafBMR","Graingrass"))

# p_variety <- ggbiplot(Strength.WR_PCA,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = Strength.WR.variety) + scale_color_manual(name = 'Varieties', values=c("blue", "orange", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
# p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



morfi.Strength.WR.section <- morfi_WRV_Sedi_na$Section



p_section <- ggbiplot(morfi_WRV_Sedi2_PCA,
                      varname.size = 3,
                      obs.scale = 0.5, 
                      ellipse=TRUE, 
                      ellipse.prob = 0.68, 
                      alpha = 0.5, 
                      groups = morfi.Strength.WR.section) + 
  scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_section + theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill="white",colour = "black", size = 2)) +
  ggtitle("PCA: MorFi + Fibre-Water Interaction Metrics")



fviz_pca_var(morfi_WRV_Sedi2_PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)




################### MorFi PCA: Low Energy ###################

morfi_WRV_Sedi_na_L <- filter(morfi_WRV_Sedi_na, grepl("L",HPH))

morfi_WRV_Sedi2_L <- subset(morfi_WRV_Sedi_na_L,select = -c(ID,Sample_name,Variety,Replicate,Section,HPH,Morfi.Rep))

morfi_WRV_Sedi2_L_PCA <- prcomp(~ ., data=morfi_WRV_Sedi2_L, na.action=na.omit, center = TRUE, scale=TRUE)




summary(morfi_WRV_Sedi2_L_PCA)

fviz_eig(morfi_WRV_Sedi2_L_PCA)

ggbiplot(morfi_WRV_Sedi2_L_PCA)



morfi.Strength.WR.section_L <- morfi_WRV_Sedi_na_L$Section



p_section <- ggbiplot(morfi_WRV_Sedi2_L_PCA,
                      varname.size = 3,
                      obs.scale = 0.5, 
                      ellipse=TRUE, 
                      ellipse.prob = 0.68, 
                      alpha = 0.5, 
                      groups = morfi.Strength.WR.section_L) + 
  scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_section + theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill="white",colour = "black", size = 2)) +
  ggtitle("PCA: MorFi + Fibre-Water Interaction Metrics")







################### MorFi PCA: Medium Energy ###################

morfi_WRV_Sedi_na_M <- filter(morfi_WRV_Sedi_na, grepl("M",HPH))

morfi_WRV_Sedi2_M <- subset(morfi_WRV_Sedi_na_M,select = -c(ID,Sample_name,Variety,Replicate,Section,HPH,Morfi.Rep))

morfi_WRV_Sedi2_M_PCA <- prcomp(~ ., data=morfi_WRV_Sedi2_M, na.action=na.omit, center = TRUE, scale=TRUE)




summary(morfi_WRV_Sedi2_M_PCA)

fviz_eig(morfi_WRV_Sedi2_M_PCA)

ggbiplot(morfi_WRV_Sedi2_M_PCA)


morfi.Strength.WR.section_M <- morfi_WRV_Sedi_na_M$Section


p_section <- ggbiplot(morfi_WRV_Sedi2_M_PCA,
                      varname.size = 3,
                      obs.scale = 0.5, 
                      ellipse=TRUE, 
                      ellipse.prob = 0.68, 
                      alpha = 0.5, 
                      groups = morfi.Strength.WR.section_M) + 
  scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_section + theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill="white",colour = "black", size = 2)) +
  ggtitle("PCA: MorFi + Fibre-Water Interaction Metrics")









################### MorFi PCA: High Energy ###################

morfi_WRV_Sedi_na_H <- filter(morfi_WRV_Sedi_na, grepl("H",HPH))

morfi_WRV_Sedi2_H <- subset(morfi_WRV_Sedi_na_H,select = -c(ID,Sample_name,Variety,Replicate,Section,HPH,Morfi.Rep))

morfi_WRV_Sedi2_H_PCA <- prcomp(~ ., data=morfi_WRV_Sedi2_H, na.action=na.omit, center = TRUE, scale=TRUE)




summary(morfi_WRV_Sedi2_H_PCA)

fviz_eig(morfi_WRV_Sedi2_H_PCA)

ggbiplot(morfi_WRV_Sedi2_H_PCA)


morfi.Strength.WR.section_H <- morfi_WRV_Sedi_na_H$Section


p_section <- ggbiplot(morfi_WRV_Sedi2_H_PCA,
                      varname.size = 3,
                      obs.scale = 0.5, 
                      ellipse=TRUE, 
                      ellipse.prob = 0.68, 
                      alpha = 0.5, 
                      groups = morfi.Strength.WR.section_H) + 
  scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_section + theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(fill="white",colour = "black", size = 2)) +
  ggtitle("PCA: MorFi + Fibre-Water Interaction Metrics")









################### Dendrogram: Morfi + Fibre-Water ###################


morfi_WRV_Sedi_means <- morfi_WRV_Sedi %>% 
  group_by(ID) %>% 
  summarise(across(.cols=where(is.numeric), mean,na.rm = TRUE))  

morfi_WRV_Sedi_means <- morfi_WRV_Sedi_means %>% dplyr:: select(1,4:26)

morfi_WRV_Sedi_means <- as.data.frame(morfi_WRV_Sedi_means)

rownames(morfi_WRV_Sedi_means) <- morfi_WRV_Sedi_means$ID
morfi_WRV_Sedi_means <- morfi_WRV_Sedi_means[,-1]

morfi_WRV_Sedi_means <- scale(morfi_WRV_Sedi_means,center=TRUE,scale=TRUE)

biomass.dist <- dist(morfi_WRV_Sedi_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"ward.D2")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

# ggdendrogram(biomass.hc.opt,rotate = TRUE)
# 
# 
# fviz_nbclust(morfi_WRV_Sedi_means, FUN = hcut, method = "wss")
# fviz_nbclust(morfi_WRV_Sedi_means, FUN = hcut, method = "silhouette")
# 
# gap_stat <- clusGap(morfi_WRV_Sedi_means, FUN = hcut, nstart = 25, K.max = 10, B = 50)
# fviz_gap_stat(gap_stat)
# 
# 
# 
# plot(biomass.hc, cex=0.6)
# rect.hclust(biomass.hc, k = 4, border = 2:5)
# 
# plot(biomass.hc.opt, cex=0.6)
# rect.hclust(biomass.hc.opt, k = 4, border = 2:5)
# 
# 
# biomass.subgroup <- cutree(biomass.hc, k = 4)
# 
# fviz_cluster(list(data = morfi_WRV_Sedi_means, cluster = biomass.subgroup),
#              ggtheme = theme_bw())




morfi.dist <- dist(t(morfi_WRV_Sedi_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"ward.D2")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)

fviz_nbclust(t(morfi_WRV_Sedi_means), FUN = hcut, method = "wss")
fviz_nbclust(t(morfi_WRV_Sedi_means), FUN = hcut, method = "silhouette")


morfi.subgroup <- cutree(morfi.hc, k = 5)

fviz_cluster(list(data = t(morfi_WRV_Sedi_means), cluster = morfi.subgroup),
             ggtheme = theme_bw())







################### Dendrogram: Nanopaper + Fibre-Water ###################

Nanopaper_WRV_Sedi_means <- Nanopaper_WRV_Sedi %>% 
  group_by(ID) %>% 
  summarise(across(.cols=where(is.numeric), mean,na.rm = TRUE))  

Nanopaper_WRV_Sedi_means <- Nanopaper_WRV_Sedi_means %>% dplyr:: select(1,3:19)

Nanopaper_WRV_Sedi_means <- as.data.frame(Nanopaper_WRV_Sedi_means)

rownames(Nanopaper_WRV_Sedi_means) <- Nanopaper_WRV_Sedi_means$ID
Nanopaper_WRV_Sedi_means <- Nanopaper_WRV_Sedi_means[,-1]

Nanopaper_WRV_Sedi_means <- scale(Nanopaper_WRV_Sedi_means,center=TRUE,scale=TRUE)

biomass.dist <- dist(Nanopaper_WRV_Sedi_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"ward.D2")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)



nanopaper.dist <- dist(t(Nanopaper_WRV_Sedi_means),"euclidean")

nanopaper.hc <- hclust(nanopaper.dist,"ward.D2")

nanopaper.hc.opt <- reorder.hclust(nanopaper.hc,nanopaper.dist)

ggdendrogram(nanopaper.hc.opt,rotate = TRUE)

fviz_nbclust(t(Nanopaper_WRV_Sedi_means), FUN = hcut, method = "wss")
fviz_nbclust(t(Nanopaper_WRV_Sedi_means), FUN = hcut, method = "silhouette")


nanopaper.subgroup <- cutree(nanopaper.hc, k = 3)

fviz_cluster(list(data = t(Nanopaper_WRV_Sedi_means), cluster = nanopaper.subgroup),
             ggtheme = theme_bw())





