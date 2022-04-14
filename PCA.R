
## PCA of MorFi data


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
library(reactablefmtr)



morfi_biomass <- read_csv("morfi_biomass.csv")

morfi_cnf<- read_csv("morfi_cnf.csv")

## Datacamp tutorial _ MorFi

morfi_cnf_PCA <- subset(morfi_cnf,select = -c(Sample_name,Replicate,Morfi.Rep,Variety,Section,HPH))


PCA.test <- prcomp(~ ., data=morfi_cnf_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.test)

fviz_eig(PCA.test)

ggbiplot(PCA.test)




################# PCA - MorFi grouped by factors ###################

  ### PCA by variety

PCA.test.variety <- c(rep("Sugargraze",48),rep("Yemen",48),rep("GreenleafBMR",36),rep("Graingrass",35))
PCA.test.variety <- factor(PCA.test.variety, levels = c("Sugargraze","Yemen","GreenleafBMR","Graingrass"))

p_variety <- ggbiplot(PCA.test,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.variety) + scale_color_manual(name = 'Varieties', values=c("blue", "orange", "darkgreen","purple")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

  ### PCA by section

PCA.test.section <- c(rep("Leaf",12),rep("Sheath",12),rep("Stem<1m",12),rep("Stem>1m",12),rep("Leaf",12),rep("Sheath",12),rep("Stem<1m",12),rep("Stem>1m",12),rep("Leaf",11),rep("Sheath",12),rep("Stem",12),rep("Leaf",12),rep("Sheath",12),rep("Stem",12))


#### Morfi PCA - Low only #########

morfi_cnf_Low <- filter(morfi_cnf, HPH =="L")
morfi_cnf_Low_PCA <- subset(morfi_cnf_Low,select = -c(Sample_name,Replicate,Morfi.Rep,Variety,Section,HPH))
PCA.test_Low <- prcomp(~ ., data=morfi_cnf_Low_PCA, na.action=na.omit, center = TRUE, scale=TRUE)
PCA.test.section_Low <- morfi_cnf_Low$Section
PCA.test.variety_Low <- morfi_cnf_Low$Variety

p_section <- ggbiplot(PCA.test_Low,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.section_Low) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_section + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

p_variety <- ggbiplot(PCA.test_Low,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.variety_Low) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

##################################################################################

#### Morfi PCA - Medium only #########

morfi_cnf_Medium <- filter(morfi_cnf, HPH =="M")
morfi_cnf_Medium_PCA <- subset(morfi_cnf_Medium,select = -c(Sample_name,Replicate,Morfi.Rep,Variety,Section,HPH))
PCA.test_Medium <- prcomp(~ ., data=morfi_cnf_Medium_PCA, na.action=na.omit, center = TRUE, scale=TRUE)
PCA.test.section_Medium <- morfi_cnf_Medium$Section
PCA.test.variety_Medium <- morfi_cnf_Medium$Variety


p_section <- ggbiplot(PCA.test_Medium,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.section_Medium) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_section + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

p_variety <- ggbiplot(PCA.test_Medium,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.variety_Medium) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

##################################################################################



#### Morfi PCA - High only #########

morfi_cnf_High <- filter(morfi_cnf, HPH =="H")
morfi_cnf_High <- morfi_cnf_High[-c(36), ]
morfi_cnf_High_PCA <- subset(morfi_cnf_High,select = -c(Sample_name,Replicate,Morfi.Rep,Variety,Section,HPH))
PCA.test_High <- prcomp(~ ., data=morfi_cnf_High_PCA, na.action=na.omit, center = TRUE, scale=TRUE)
PCA.test.section_High <- morfi_cnf_High$Section
PCA.test.variety_High <- morfi_cnf_High$Variety


p_section <- ggbiplot(PCA.test_High,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.section_High) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_section + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

p_variety <- ggbiplot(PCA.test_High,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.variety_High) + scale_color_discrete(name = 'Sections') +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

##################################################################################





p_section <- ggbiplot(PCA.test,varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.section) + scale_color_discrete(name = 'Sections') + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_section + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

  ### PCA by HPH

PCA.test.HPH <- c(rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",3),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4),rep("L",4),rep("M",4),rep("H",4))
PCA.test.HPH <- factor(PCA.test.HPH, levels = c("L","M","H"))

#PCAs 1 & 2
p_HPH12 <- ggbiplot(PCA.test,choices=1:2, varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.HPH) + scale_color_manual(name = 'HPH Energy Level', values=c("lightgreen", "orange", "red")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_HPH12 + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

#PCAs 2 & 3
p_HPH23 <- ggbiplot(PCA.test,choices=2:3, varname.size = 3,obs.scale = 0.5, ellipse=TRUE, ellipse.prob = 0.68, alpha = 0.5, groups = PCA.test.HPH) + scale_color_manual(name = 'HPH Energy Level', values=c("lightgreen", "orange", "red")) + theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank())+geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
p_HPH23 + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))





################# PCA - Nanopaper grouped by factors ###################

# Datacamp tutorial _ Nanopaper_all metrics

nanopaper_PCA <- subset(nanopaper,select = -c(biomass_ID,sample_name,rep,variety,section,HPH,nanopaper_rep,strain,UTS,modulus,toughness,TI))
colnames(nanopaper_PCA)[6] <- "strain"
colnames(nanopaper_PCA)[7] <- "UTS"
colnames(nanopaper_PCA)[8] <- "modulus"
colnames(nanopaper_PCA)[9] <- "toughness"
colnames(nanopaper_PCA)[10] <- "TI"

# write.csv(x=nanopaper_PCA, file='strength.v.WRV_PCA.csv')



PCA.nanopaper <- prcomp(~ ., data=nanopaper_PCA, center = TRUE, scale=TRUE)
summary(PCA.nanopaper)
a <- summary(PCA.nanopaper)$importance 

write.csv(x=a, file='a.csv')

fviz_eig(PCA.nanopaper)

ggbiplot(PCA.nanopaper)

  ### PCA by HPH

PCA.nanopaper.HPH <- c(rep("L",7),rep("M",8),rep("H",12),rep("L",16),rep("M",15),rep("H",16),rep("L",16),rep("M",16),rep("H",15),rep("L",15),rep("M",15),rep("H",15),
                       rep("L",13),rep("M",16),rep("H",11),rep("L",16),rep("M",16),rep("H",15),rep("L",16),rep("M",16),rep("H",16),rep("L",15),rep("M",14),rep("H",15),
                       rep("L",8),rep("M",12),rep("H",11),rep("L",13),rep("M",13),rep("H",12),rep("L",14),rep("M",16),rep("H",11),
                       rep("L",14),rep("M",14),rep("H",14),rep("L",12),rep("M",15),rep("H",14),rep("L",15),rep("M",14),rep("H",15))

PCA.nanopaper.HPH <- factor(PCA.nanopaper.HPH, levels = c("L","M","H"))
PCA.nanopaper.reverse <- PCA.nanopaper
PCA.nanopaper.reverse$x[,1] <- - PCA.nanopaper.reverse$x[,1]
PCA.nanopaper.reverse$rotation[,1] <- - PCA.nanopaper.reverse$rotation[,1]


p_nanopaper_HPH <- ggbiplot(PCA.nanopaper.reverse,varname.size = 4,ellipse=TRUE, obs.scale = 0.5, alpha = 0.5, groups = PCA.nanopaper.HPH) + 
  scale_color_manual(name = 'HPH Energy Level', values=c("lightgreen", "orange", "red")) + 
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) +
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") +
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_nanopaper_HPH + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))


  ### PCA by Section

PCA.nanopaper.section <- c(rep("Leaf",27),rep("Sheath",47),rep("Stem<1m",47),rep("Stem>1m",45),
                           rep("Leaf",40),rep("Sheath",47),rep("Stem<1m",48),rep("Stem>1m",44),
                           rep("Leaf",31),rep("Sheath",38),rep("Stem",41),
                           rep("Leaf",42),rep("Sheath",41),rep("Stem",44))

p_nanopaper_section <- ggbiplot(PCA.nanopaper,varname.size = 4,ellipse=TRUE, obs.scale = 0.5, alpha = 0.5, groups = PCA.nanopaper.section)+
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) +
  scale_color_discrete(name = 'Sections') +
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") +
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_nanopaper_section + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))

  ### PCA by Variety

PCA.nanopaper.variety <- c(rep("Sugargraze",166),
                           rep("Yemen",179),
                           rep("GreenleafBMR",110),
                           rep("Graingrass",127))

p_nanopaper_variety <- ggbiplot(PCA.nanopaper,varname.size = 4,ellipse=TRUE, obs.scale = 0.5, alpha = 0.4,groups = PCA.nanopaper.variety) +
  theme(legend.direction = 'horizontal', legend.position = 'top', panel.background = element_blank()) +
  scale_color_manual(name = 'Varieties', values=c("blue", "orange", "darkgreen","purple")) +
  geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed") +
  geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")

p_nanopaper_variety + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill="white",colour = "black", size = 2))



## Datacamp tutorial _ Nanopaper_4 metrics


NP_PCA <- subset(nanopaper,select = -c(ID,sample_name,rep,variety,section,HPH,nanopaper_rep,weight,thickness,density,porosity,grammage,strain,UTS,modulus,toughness,TI))
PCA.NP <- prcomp(~ ., data=NP_PCA, center = TRUE, scale=TRUE)

# summary(PCA.NP)
# 
# fviz_eig(PCA.NP)
# 
# ggbiplot(PCA.NP)
# 
# 
# PCA.NP.HPH <- c(rep("L",7),rep("M",8),rep("H",12),rep("L",16),rep("M",15),rep("H",16),rep("L",16),rep("M",16),rep("H",15),rep("L",15),rep("M",15),rep("H",15),
#                        rep("L",13),rep("M",16),rep("H",11),rep("L",16),rep("M",16),rep("H",15),rep("L",16),rep("M",16),rep("H",16),rep("L",15),rep("M",14),rep("H",15),
#                        rep("L",8),rep("M",12),rep("H",11),rep("L",13),rep("M",13),rep("H",12),rep("L",14),rep("M",16),rep("H",11),
#                        rep("L",14),rep("M",14),rep("H",14),rep("L",12),rep("M",15),rep("H",14),rep("L",15),rep("M",14),rep("H",15))
# ggbiplot(PCA.NP,ellipse=TRUE, groups = PCA.NP.HPH)
# 
# 
# PCA.NP.section <- c(rep("Leaf",27),rep("Sheath",47),rep("Stem<1m",47),rep("Stem>1m",45),
#                            rep("Leaf",40),rep("Sheath",47),rep("Stem<1m",48),rep("Stem>1m",44),
#                            rep("Leaf",31),rep("Sheath",38),rep("Stem",41),
#                            rep("Leaf",42),rep("Sheath",41),rep("Stem",44))
# ggbiplot(PCA.NP,ellipse=TRUE, groups = PCA.NP.section)
# 
# 
# PCA.NP.variety <- c(rep("Sugargraze",166),
#                            rep("Yemen",179),
#                            rep("GreenleafBMR",110),
#                            rep("Graingrass",127))
# ggbiplot(PCA.NP,ellipse=TRUE, groups = PCA.NP.variety)




################# PCA - Extrusion ###################


## Datacamp tutorial _ Extrusion.NP

Ex.NP_PCA <- subset(morfi_ex,select = c(strain_star,E_star,UT_star,TI_star))
PCA.Ex.NP <- prcomp(~ ., data=Ex.NP_PCA, center = TRUE, scale=TRUE)
summary(PCA.Ex.NP)

ggbiplot(PCA.Ex.NP)


PCA.Ex.NP.pass <- c(rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4))
ggbiplot(PCA.Ex.NP,ellipse=TRUE, groups = PCA.Ex.NP.pass)


PCA.Ex.NP.SC <- c(rep("10",12),rep("15",12),rep("20",12),rep("25",12))
ggbiplot(PCA.Ex.NP,ellipse=TRUE, groups = PCA.Ex.NP.SC)



## Datacamp tutorial _ Extrusion.Morfi

morfi_ex_PCA <- subset(morfi_ex,select = -c(sample_name,SC,pass,strain_star,E_star,UTS_star,UT_star,TI_star))

PCA.Ex.Morfi <- prcomp(~ ., data=morfi_ex_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.Ex.Morfi)

ggbiplot(PCA.Ex.Morfi, loadings.label.repel=TRUE, obs.scale = 0.8, var.scale = 0.8)


PCA.Ex.Morfi.pass <- c(rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4))
ggbiplot(PCA.Ex.Morfi,ellipse=TRUE, groups = PCA.Ex.Morfi.pass)


PCA.Ex.Morfi.SC <- c(rep("10",12),rep("15",12),rep("20",12),rep("25",12))
ggbiplot(PCA.Ex.Morfi,ellipse=TRUE, groups = PCA.Ex.Morfi.SC,obs.scale = 0.5)




## Datacamp tutorial _ Extrusion.Morfi.Old vs New

morfi_ex_oldnew_PCA <- subset(morfi_ex_oldnew,select = -c(sample_name,age,SC,pass))

PCA.Morfi_Ex_oldnew <- prcomp(~ ., data=morfi_ex_oldnew_PCA, na.action=na.omit, center = TRUE, scale=TRUE)

summary(PCA.Morfi_Ex_oldnew)

ggbiplot(PCA.Morfi_Ex_oldnew, obs.scale = 0.8, var.scale = 0.8)

morfi_ex_oldnew_PCA.age <- c(rep("new",28),rep("old",28))
ggbiplot(PCA.Morfi_Ex_oldnew, ellipse=TRUE, obs.scale = 1.5, var.scale = 0.5, groups = morfi_ex_oldnew_PCA.age)


p_oldnew <- ggbiplot(PCA.Morfi_Ex_oldnew,choices=1:2, 
                     varname.size = 3,
                     obs.scale = 1, 
                     ellipse=TRUE, 
                     ellipse.prob = 0.68, 
                     alpha = 0.5, 
                     groups = morfi_ex_oldnew_PCA.age) + 
  scale_color_manual(name = 'Pulp age', values=c("#333333", "orange")) + 
  theme(legend.direction = 'horizontal', 
        legend.position = 'top', 
        panel.background = element_blank()) +
  geom_hline(yintercept=0,
             colour = "grey80", 
             size = 0.5, 
             lty="dashed") + 
  geom_vline(xintercept=0,
             colour = "grey80", 
             size = 0.5, 
             lty="dashed")

p_oldnew + theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.background = element_rect(fill="white",colour = "black", size = 2))





morfi_ex_oldnew_PCA.pass <- c(rep("1",4),rep("2",4),rep("3",4),rep("5",4),rep("1",4),rep("2",4),rep("3",4),rep("1",4),rep("2",4),rep("3",4),rep("5",4),rep("1",4),rep("2",4),rep("3",4))
ggbiplot(PCA.Morfi_Ex_oldnew,ellipse=TRUE, groups = morfi_ex_oldnew_PCA.pass)


morfi_ex_oldnew_PCA.SC <- c(rep("20",16),rep("25",12),rep("20",16),rep("25",12))
ggbiplot(PCA.Morfi_Ex_oldnew,ellipse=TRUE, groups = morfi_ex_oldnew_PCA.SC)




# ## CRAN tutorial
#
# fviz_eig(PCA.test)
#
# fviz_pca_var(PCA.test,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )



    #Dendrogram TEST

# morfi_cnf_dend <- subset(morfi_cnf,select = -c(Sample_name,Replicate,Morfi.Rep))
# morfi_cnf_dend$name <- str_c(morfi_cnf_dend$Variety,"-",morfi_cnf_dend$Section,"-",morfi_cnf_dend$HPH)
# 
# t_morfi_cnf_dend <- as.data.frame(as.matrix(t(morfi_cnf_PCA[,1:21])))
# 
# colnames(t_morfi_cnf_dend) <- morfi_cnf_dend$name
# 
# 
# #######Trying to get the sample names as rows: rownames(morfi_cnf_PCA) <- morfi_cnf_dend$name
# 
# # rownames(t_morfi_cnf_dend) <- colnames(morfi_cnf_PCA)
# #
# # rownames(morfi_cnf_PCA) <- colnames(morfi_cnf_PCA)
# 
# 
# 
# res.pca <- PCA(t_morfi_cnf_dend[,-108], ncp = 4, graph = TRUE,na.rm=TRUE)
# 
# res.pca <- prcomp(~ ., data=t_morfi_cnf_dend, na.action=na.omit, center = TRUE, scale=TRUE)
# 
# res.hcpc <- HCPC(res.pca, graph = TRUE)
# 
# fviz_dend(res.hcpc,
#           cex = 0.8,                     # Label size
#           palette = "jco",               # Color palette see ?ggpubr::ggpar
#           rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#           rect_border = "jco",           # Rectangle color
#           labels_track_height = 0,       # Augment the room for labels
#           main = "Dendrogram",
#           xlab = "Morfi Parameters", ylab = "Distance", sub = "...",
#           horiz = TRUE,
#           repel = TRUE,
#           show_labels = TRUE,
#           color_labels_by_k = TRUE)





##################### Dendrogram - MorFi #####################


  #### Biomass sample Dendrogram

morfi_cnf$biomass_ID <- with(morfi_cnf,paste(Variety,Section,HPH, sep= "-"))

morfi_cnf_means <- morfi_cnf %>% group_by(biomass_ID) %>% summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
morfi_cnf_means <- morfi_cnf_means %>% dplyr:: select(1,4:24)
  
morfi_cnf_means <- as.data.frame(morfi_cnf_means)

rownames(morfi_cnf_means) <- morfi_cnf_means$biomass_ID
morfi_cnf_means <- morfi_cnf_means[,-1]

morfi_cnf_means <- scale(morfi_cnf_means,center=TRUE,scale=TRUE)

biomass.dist <- dist(morfi_cnf_means,"euclidean")

biomass.hc <- hclust(biomass.dist,"ward.D2")

biomass.hc.opt <- reorder.hclust(biomass.hc,biomass.dist)

ggdendrogram(biomass.hc.opt,rotate = TRUE)


fviz_nbclust(morfi_cnf_means, FUN = hcut, method = "wss")
fviz_nbclust(morfi_cnf_means, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(morfi_cnf_means, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)



plot(biomass.hc, cex=0.6)
rect.hclust(biomass.hc, k = 4, border = 2:5)

plot(biomass.hc.opt, cex=0.6)
rect.hclust(biomass.hc.opt, k = 4, border = 2:5)


biomass.subgroup <- cutree(biomass.hc, k = 4)

fviz_cluster(list(data = morfi_cnf_means, cluster = biomass.subgroup),
             ggtheme = theme_bw())



# https://www.datanovia.com/en/blog/k-means-clustering-visualization-in-r-step-by-step-guide/
# set.seed(123)
# res.km <- kmeans(scale(morfi_cnf_PCA), 3, nstart = 25)
# 
# ind.coord <- as.data.frame(get_pca_ind(PCA.test)$coord)
# ind.coord$cluster <- factor(morfi_cnf_means$cluster)





# morfi_cnf_PCA_means <- morfi_cnf_PCA %>% 
#   dplyr::group_by(Variety,Section,HPH) %>% 
#   dplyr::summarise(across(.cols=everything(),mean,na.rm=TRUE)) %>% 
#   dplyr::ungroup()
# 
# df <- morfi_cnf_PCA
# df <- na.omit(df)
# df <- scale(df)
# morfi_cnf <- na.omit(morfi_cnf)
# rownames(df) <- morfi_cnf$biomass_ID



  #### Morfi parameters Dendrogram

morfi.dist <- dist(t(morfi_cnf_means),"euclidean")

morfi.hc <- hclust(morfi.dist,"ward.D2")

morfi.hc.opt <- reorder.hclust(morfi.hc,morfi.dist)

ggdendrogram(morfi.hc.opt,rotate = TRUE)

fviz_nbclust(t(morfi_cnf_means), FUN = hcut, method = "wss")
fviz_nbclust(t(morfi_cnf_means), FUN = hcut, method = "silhouette")


morfi.subgroup <- cutree(morfi.hc, k = 5)

fviz_cluster(list(data = t(morfi_cnf_means), cluster = morfi.subgroup),
             ggtheme = theme_bw())





##################### Dendrogram - Nanopaper #####################


#### Biomass sample Dendrogram

nanopaper$biomass_ID <- with(nanopaper,paste(variety,section,HPH, sep= "-"))

nanopaper_means <- nanopaper %>% group_by(biomass_ID) %>% summarise(across(.cols=where(is.numeric),mean,na.rm = TRUE))  
nanopaper_means <- nanopaper_means %>% dplyr:: select(1,5:19)

nanopaper_means <- as.data.frame(nanopaper_means)

rownames(nanopaper_means) <- nanopaper_means$biomass_ID
nanopaper_means <- nanopaper_means[,-1]

nanopaper_means <- scale(nanopaper_means,center=TRUE,scale=TRUE)

nanopaper.biomass.dist <- dist(nanopaper_means,"euclidean")

nanopaper.biomass.hc <- hclust(nanopaper.biomass.dist,"ward.D2")

nanopaper.biomass.hc.opt <- reorder.hclust(nanopaper.biomass.hc,nanopaper.biomass.dist)

plot(nanopaper.biomass.hc, cex=0.6)
rect.hclust(nanopaper.biomass.hc, k = 5, border = 2:5)

ggdendrogram(nanopaper.biomass.hc.opt,rotate = TRUE)


fviz_nbclust(nanopaper_means, FUN = hcut, method = "wss")
fviz_nbclust(nanopaper_means, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(nanopaper_means, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

nanopaper.biomass.subgroup <- cutree(nanopaper.biomass.hc, k = 4)

fviz_cluster(list(data = nanopaper_means, cluster = nanopaper.biomass.subgroup),
             ggtheme = theme_bw())



#### Nanopaper metric Dendrogram

nanopaper.metric.dist <- dist(t(nanopaper_means),"euclidean")

nanopaper.metric.hc <- hclust(nanopaper.metric.dist,"ward.D2")

nanopaper.metric.hc.opt <- reorder.hclust(nanopaper.metric.hc,nanopaper.metric.dist)

ggdendrogram(nanopaper.metric.hc.opt,rotate = TRUE)

fviz_nbclust(t(nanopaper_means), FUN = hcut, method = "wss")
fviz_nbclust(t(nanopaper_means), FUN = hcut, method = "silhouette")

gap_stat <- clusGap(t(nanopaper_means), FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


nanopaper.metric.subgroup <- cutree(nanopaper.metric.hc, k = 4)

fviz_cluster(list(data = t(nanopaper_means), cluster = nanopaper.metric.subgroup),
             ggtheme = theme_bw())




##################### Tanglegram #####################


nanopaper.biomass.hc <- hclust(nanopaper.biomass.dist,"ward.D2")
morfi.biomass.hc <- hclust(biomass.dist,"ward.D2")

plot(nanopaper.biomass.hc, cex = 0.8)
plot(morfi.biomass.hc, cex = 0.8)

dend_nanopaper <- as.dendrogram (nanopaper.biomass.hc)
dend_morfi <- as.dendrogram (morfi.biomass.hc)

tanglegram(dend_nanopaper, dend_morfi)

dend_list <- dendlist(dend_nanopaper, dend_morfi)

# all.equal(dend_nanopaper, dend_morfi)

tanglegram(dend_nanopaper, dend_morfi, sort = TRUE,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = TRUE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)






##################### Customised Tanglegram #####################

dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}

hcdata1 <- dendro_data_k(nanopaper.biomass.hc, 5)
hcdata2 <- dendro_data_k(morfi.biomass.hc, 5)

cols    <- c("#a9a9a9", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

p1 <- plot_ggdendro(hcdata1,
                    direction   = "lr",
                    scale.color = cols,
                    expand.y    = 0.2) +
  theme_void()

p2 <- plot_ggdendro(hcdata2,
                    direction   = "rl",
                    scale.color = cols,
                    expand.y    = 0.2) +
  theme_void()

idx <- data.frame(y1 = 1:nrow(hcdata1$labels),
                  y2 = match(hcdata1$labels$label, hcdata2$labels$label))

p3 <- ggplot() +
  geom_segment(data     = idx, 
               aes(x    = 0,
                   y    = y1,
                   xend = 1,
                   yend = y2),
               color    = "grey") +
  theme_void()

grid.arrange(p1, p3, p2, ncol = 3, widths = c(2, 1, 2),
             top = textGrob(paste("Entanglement =", round(entanglement(dend_list), 2)),
                            gp=gpar(fontsize=20,font=3)))



##################### Tanglegram Correlations #####################


cor_bakers_gamma(dend_nanopaper, dend_morfi)



# # The function cor.dendlist() is used to compute “Baker” or “Cophenetic”
# correlation matrix between a list of trees. The value can range between -1 to 1.
# With near 0 values meaning that the two trees are not statistically similar.

cor.dendlist(dend_list, method = "cophenetic")
cor_cophenetic(dend_nanopaper, dend_morfi)
cor.dendlist(dend_list, method = "baker")
cor_bakers_gamma(dend_nanopaper, dend_morfi)


dend_list_names <- dendlist("Nanopaper" = dend_nanopaper, "MorFi" = dend_morfi)
cors <- cor.dendlist(dend_list_names)
round(cors, 2)

# corrplot(cors, "pie", "lower")










## Significance Testing w/ PCA for CNF pulp

##assign PC cutoff value to select number of components/variables included


parameters <- as.data.frame(abs(PCA.test$rotation[,1:6]))

write.csv(x=parameters, file='parameters.csv')

# highlight_min_max(
#   parameters,
#   min_font_color = "red",
#   max_font_color = "green",
#   min_highlighter = NULL,
#   max_highlighter = NULL
# )




v1 <- max(abs(PCA.test$rotation[,1]),na.rm=TRUE)
v2 <- max(abs(PCA.test$rotation[,2]),na.rm=TRUE)
v3 <- max(abs(PCA.test$rotation[,3]),na.rm=TRUE)
v4 <- max(abs(PCA.test$rotation[,4]),na.rm=TRUE)
# v5 <- max(abs(PCA.test$rotation[,5]),na.rm=TRUE)
# v6 <- max(abs(PCA.test$rotation[,6]),na.rm=TRUE)
# v7 <- max(abs(PCA.test$rotation[,7]),na.rm=TRUE)

var1 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,1])==v1]
var2 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,2])==v2]
var3 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,3])==v3]
var4 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,4])==v4]
# var5 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,5])==v5]
# var6 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,6])==v6]
# var7 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,7])==v7]


    #Variable 1 CLD plot

aov.var1 <- aov(as.formula(paste(eval(var1),"HPH*Variety/Section",sep="~")), data = morfi_cnf)

summary(aov.var1)

MM.aov.var1 <- emmeans(aov.var1,specs=~HPH:Variety:Section)

cld.var1 <- cld(MM.aov.var1,reversed=TRUE,Letters=letters)

cld.var1$ID <- str_c(cld.var1$Variety,'-',cld.var1$Section,'-',cld.var1$HPH)
cld.var1$ID <- as.factor(cld.var1$ID)

cld.var1 %>%
  mutate(ID = fct_reorder(ID, emmean,.fun='mean')) %>%
  ggplot(aes(x=emmean, y=ID))+
  geom_errorbarh(aes(xmin=lower.CL,
                    xmax=upper.CL),
                size=3,
                colour = 'blue',
                alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  geom_text(aes(x=30,label=cld.var1$.group),
            fontface=2,
            hjust = 1) +
  xlab("Mean Fibre Width (um)") +
  ylab("Sample ID")










    ##Testing for one section

aov.var1.section <- aov(as.formula(paste(eval(var1),"type*Variety",sep="~")), data = morfi_cnf,subset=Section=="Stem(<1m)")
summary(aov.var1.section)
MM.aov.var1.section <- emmeans(aov.var1.section,specs=~type:Variety)
cld.stem_1m <- cld(MM.aov.var1.section,reversed=TRUE,Letters=letters)


# cld.stem_1m %>%
#   mutate(ID = fct_reorder(ID, emmean,.fun='mean')) %>%
#   ggplot(aes(x=emmean, y=ID))+
#   geom_errorbarh(aes(xmin=lower.CL,
#                      xmax=upper.CL),
#                  size=3,
#                  colour = 'blue',
#                  alpha=0.3)+
#   geom_point(stat="identity", size=2.5)+
#   geom_text(aes(x=30,label=cld.var1$.group),
#             fontface=2,
#             hjust = 1) +
#   xlab("Mean Fibre Width (um)") +
#   ylab("Sample ID")








# 
# ############ Significance Ranking w/ PCA - Biomass - Fibre Width ############
# 
# morfi_biomass <- read_csv("morfi_biomass.csv")
# colnames(morfi_biomass)[1] <- c('sample_name')
# morfi_biomass_PCA <- subset(morfi_biomass,select = -c(sample_name,Replicate,Morfi.Rep,Variety,Section))
# 
# 
# PCA.biomass <- prcomp(~ ., data=morfi_biomass_PCA, na.action=na.omit, center = TRUE, scale=TRUE)
# 
# v1 <- max(abs(PCA.biomass$rotation[,1]),na.rm=TRUE)
# v2 <- max(abs(PCA.biomass$rotation[,2]),na.rm=TRUE)
# v3 <- max(abs(PCA.biomass$rotation[,3]),na.rm=TRUE)
# v4 <- max(abs(PCA.biomass$rotation[,4]),na.rm=TRUE)
# 
# var1 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,1])==v1]
# var2 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,2])==v2]
# var3 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,3])==v3]
# var4 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,4])==v4]
# 
# biomass.aov.var1 <- aov(as.formula(paste(eval(var1),"Variety/Section",sep="~")), data = morfi_biomass)
# 
# MM.biomass.aov.var1 <- emmeans(biomass.aov.var1,specs=~Variety:Section)
# 
# cld.biomass.var1 <- cld(MM.biomass.aov.var1,reversed=TRUE,Letters=letters)
# 
# cld.biomass.var1$ID <- str_c(cld.biomass.var1$Variety,'-',cld.biomass.var1$Section)
# cld.biomass.var1$ID <- as.factor(cld.biomass.var1$ID)
# 
# cld.biomass.var1 %>%
#   mutate(ID = fct_reorder(ID, emmean,.fun='mean')) %>%
#   ggplot(aes(x=emmean, y=ID))+
#   geom_errorbar(aes(xmin=lower.CL,
#                     xmax=upper.CL),
#                 size=5,
#                 width=0,
#                 colour = 'blue',
#                 alpha=0.3)+
#   geom_point(stat="identity", size=2.5)+
#   geom_text(aes(x=55,label=cld.biomass.var1$.group),
#             fontface=2,
#             hjust = 0.5) +
#   xlab("Mean Fibre Width (um)") +
#   ylab("Sample ID")











# ############ Significance Ranking w/ PCA - Biomass - Fibre Length ############
# 
# morfi_biomass <- read_csv("morfi_biomass.csv")
# colnames(morfi_biomass)[1] <- c('sample_name')
# morfi_biomass_PCA <- subset(morfi_biomass,select = -c(sample_name,Replicate,Morfi.Rep,Variety,Section))
# 
# 
# PCA.biomass <- prcomp(~ ., data=morfi_biomass_PCA, na.action=na.omit, center = TRUE, scale=TRUE)
# 
# v1 <- max(abs(PCA.biomass$rotation[,1]),na.rm=TRUE)
# v2 <- max(abs(PCA.biomass$rotation[,2]),na.rm=TRUE)
# v3 <- max(abs(PCA.biomass$rotation[,3]),na.rm=TRUE)
# v4 <- max(abs(PCA.biomass$rotation[,4]),na.rm=TRUE)
# 
# var1 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,1])==v1]
# var2 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,2])==v2]
# var3 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,3])==v3]
# var4 <- rownames(PCA.biomass$rotation)[abs(PCA.biomass$rotation[,4])==v4]
# 
# biomass.aov.var1 <- aov(as.formula(paste(eval(var1),"Variety/Section",sep="~")), data = morfi_biomass)
# 
# MM.biomass.aov.var1 <- emmeans(biomass.aov.var1,specs=~Variety:Section)
# 
# cld.biomass.var1 <- cld(MM.biomass.aov.var1,reversed=TRUE,Letters=letters)
# 
# cld.biomass.var1$ID <- str_c(cld.biomass.var1$Variety,'-',cld.biomass.var1$Section)
# cld.biomass.var1$ID <- as.factor(cld.biomass.var1$ID)
# 
# cld.biomass.var1 %>%
#   mutate(ID = fct_reorder(ID, emmean,.fun='mean')) %>%
#   ggplot(aes(x=emmean, y=ID))+
#   geom_errorbar(aes(xmin=lower.CL,
#                     xmax=upper.CL),
#                 size=5,
#                 width=0,
#                 colour = 'blue',
#                 alpha=0.3)+
#   geom_point(stat="identity", size=2.5)+
#   geom_text(aes(x=55,label=cld.biomass.var1$.group),
#             fontface=2,
#             hjust = 0.5) +
#   xlab("Mean Fibre Width (um)") +
#   ylab("Sample ID")






















# ## Creating PLS model to predict nanopaper properties from fibre morphology data
# 
# colnames(morfi_cnf)[1] <- "Sample_name"
# 
# pls_morfi <- morfi_cnf %>% select(-c("Sample_name", "rep", "Morfi.Rep.CNF", "type"))






