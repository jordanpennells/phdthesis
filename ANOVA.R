
# Install libraries

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(effects)
library(lme4)
library(Rcpp)
library(multcomp)
library(stringr)
library(multcomp)
library(multcompView)
library(car)
library(emmeans)
library(magrittr)
library(effects)
library(rstatix)
library(report)
library(haven) # to load the SPSS .sav file
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)



rm(model)
rm(mod1)
rm(model1)


# Read CSV files

nanopaper_stem <- read_csv("nanopaper_ANOVA.csv")
nanopaper_samples <- read_csv("nanopaper_samples.csv")


nanopaper_stem_stature$HPH.Level <- factor(nanopaper_stem_stature$HPH, levels = c("L", "M", "H"))
nanopaper_stem_stature$ID = str_c(nanopaper_stem_stature$variety,'-',nanopaper_stem_stature$section,'-',nanopaper_stem_stature$HPH,'-',nanopaper_stem_stature$nanopaper,'-',nanopaper_stem_stature$strip)



# Make HPH energy an ordered factor

HPH.fac <- factor(nanopaper_samples$HPH, levels = c("L", "M", "H"))

nanopaper_stem$nanopaper <- as.factor(nanopaper_stem$nanopaper)

nanopaper_stem$strip <- as.factor(nanopaper_stem$strip)





#################### Mixed-effects Linear Model (lmer) ####################


nanopaper_stem$HPH.fac <- factor(nanopaper_stem$HPH, levels = c("L", "M", "H"))

model <-  lmer(TI_star~(variety*HPH + (variety:section*HPH) + (1|nanopaper/strip)),data=nanopaper_stem,na.action="na.omit")

# model <-  lmer(strength ~ (variety*energy + (variety:section*energy) + (1|nanopaper/strip)), data=nanopaper_samples)

summary(model)

report(model)


### ANOVA

anova(model)

### ANODEV

Anova(model,type="III",test.statistic="F")




#################### Random Model ####################

random.model <-  lmer(TI_star~(1|variety) * (1|HPH) * (1|section) + (1|nanopaper/strip),data=nanopaper_stem,na.action="na.omit")

# random.model <-  lmer(strength ~ (1|variety) * (1|HPH) * (1|section) + (1|nanopaper/strip), data = nanopaper_samples)

summary(random.model)





random.model2 <-  lmer(TI_star~(1|variety) * (1|HPH) * (1|section:variety) + (1|nanopaper/strip),data=nanopaper_stem,na.action="na.omit")

summary(random.model2)








############# Repeated measures ANOVA #############


repeated.ANOVA <- anova_test(data = nanopaper_stem_stature, 
                             formula = TI_star~(variety*HPH + (variety:section*HPH) + (1|nanopaper/strip)),
                             dv = "TI_star",
                             wid = ID,
                             within = c(Stature, section, HPH.Level))

get_anova_table(repeated.ANOVA)



############# AllEffects Graph #############

plot(allEffects(model), 
     main = "All Effects Plot for Biomass-to-Nanopaper Model", 
     xlab="Sorghum Variety", 
     ylab="Nanopaper Tensile Index (Nm/g)")









ranef(model)

# e <- allEffects(model)
# e
# plot(e)

nanopaper_stem$HPH.num <- gsub("L",0,gsub("M",0.5,gsub("H",1,nanopaper_stem$HPH)))


# ggplot(data = nanopaper_stem, 
#        aes(x   = HPH.num,
#            y   = TI_star, 
#            col = variety.fac))+
#   geom_point(size     = 1, 
#              alpha    = .7, 
#              position = "jitter")+
#   geom_smooth(method   = lm,
#               se       = T, 
#               size     = 1.5, 
#               linetype = 1, 
#               alpha    = .7)+
#   theme_bw()








##NEED TO FIX
plot(e,multiline=TRUE,confint=TRUE,ci.style="bars"
     ,main="Effect of Different Sorghum Varieties & Section on Nanopaper Strength"
     ,xlab="Sorghum Variety"
     ,ylab="Tensile Index (Nm/g)")




ggplot(nanopaper_stem,aes(x=section,y=TI_star,col=variety)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~variety) +
  theme_bw()


ggplot(nanopaper_stem,aes(x=variety,y=TI_star,col=variety)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~nanopaper) +
  theme_bw()

ggplot(nanopaper_stem,aes(x=variety,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


### Include in Thesis -> Nanopaper duplicate effect ###
ggplot(nanopaper_stem,aes(x=HPH.fac,y=TI_star,col=nanopaper)) + 
  geom_jitter(alpha=0.5) + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~variety) +
  xlab("Homogenisation Energy Level") +
  ylab("Nanopaper Tensile Index (Nm/g)") +
  theme_bw() +
  theme(legend.key.size = unit(2, 'cm'),
        legend.title = element_text(face="bold",size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14))


ggplot(nanopaper_stem,aes(x=section,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


ggplot(nanopaper_stem,aes(x=nanopaper,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


### Include in Thesis -> Nanopaper strip replicate effect ###
ggplot(nanopaper_stem,aes(x=strip,y=TI_star,col=variety)) + 
  geom_point(alpha=0.5) + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~variety) +
  xlab("Nanopaper Strip") +
  ylab("Nanopaper Tensile Index (Nm/g)") +
  theme_bw() +
  theme(legend.key.size = unit(2, 'cm'),
        legend.title = element_text(face="bold",size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14))


ggplot(nanopaper_stem,aes(x=strip,y=TI_star,col=section)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~section) +
  theme_bw()


ggplot(nanopaper_stem,aes(x=strip,y=TI_star,col=HPH)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()



ggplot(nanopaper_stem,aes(x=strip,y=TI_star,col=HPH)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


##### Add stature column #####

# nanopaper_stem_stature <- transform(nanopaper_stem, c= ifelse(variety=="Sugargraze", 
#                                                               "Tall", 
#                                                               "Short"))

myFunction <- function(x){
  variety <- x[1]
  #further values ignored (if there are more than 2 columns)
  value <- if(variety == "Sugargraze" || variety ==  "Yemen") "Tall" else "Short"
  #or more complicated stuff
  return(value)
}

nanopaper_stem_stature$Stature <- apply(nanopaper_stem_stature, 1, myFunction)


ggplot(nanopaper_stem_stature,aes(x=HPH.fac,y=TI_star,col=Stature)) + 
  geom_jitter() +
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~section) +
  theme_bw()


ggplot(nanopaper_stem_stature,aes(x=HPH.fac,y=TI_star,col=Stature)) + 
  geom_jitter() +
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~Stature) +
  theme_bw()


ggplot(nanopaper_stem_stature,aes(x=HPH.fac,y=TI_star,col=Stature)) + 
  geom_jitter() +
  geom_boxplot(alpha=0.2) + 
  facet_wrap(.~Stature + section) +
  theme_bw()

# ggplot(nanopaper_stem_stature,aes(x=HPH.fac,y=TI_star,col=Stature)) + 
#   geom_jitter() +
#   geom_boxplot(alpha=0.2) + 
#   facet_wrap(~ section, nrow =1) +
#   theme_bw()












####### Check for assumptions -> Normality & Equality of variance ##########

### Check normality assumption

ggqqplot(residuals(model), title = "QQ-plot test of normality")

ggqqplot(nanopaper_stem, "TI_star", ggtheme = theme_bw()) +
  facet_grid(variety ~ section)


shapiro.test(residuals(model))

n <- nanopaper_stem %>%
  group_by(variety, section) %>%
  shapiro_test(TI_star)

n

write.csv(x=n, file='shapiro_wilk.csv')


ggqqplot(nanopaper_stem_stature, "TI_star", ggtheme = theme_bw()) +
  facet_grid(Stature + section ~ HPH.Level, labeller = "label_both")

# ggplot(n, aes(x=))


### Homogeneity of variance assumption


plot(fitted(model), resid(model, type = "pearson"))# this will create the plot
abline(0,0, col="red")

# dplyr::


levene <- nanopaper_stem %>%
  group_by(variety) %>%
  levene_test(TI_star ~ section)

write.csv(x=levene, file='levene.csv')


nanopaper_stem %>%
  group_by(variety) %>%
  levene_test(TI_star ~ HPH)



nanopaper_stem %>%
  group_by(HPH) %>%
  levene_test(TI_star ~ section)

nanopaper_stem %>%
  group_by(HPH) %>%
  levene_test(TI_star ~ variety)

levene_test(model)

#The Levene’s test is significant (p < 0.05), therefore we cannot assume the homogeneity of variances in the different groups.
#The assumption of homogeneity of variance means that the level of variance for a particular variable is constant across the sample.
# In ANOVA, when homogeneity of variance is violated there is a greater probability of falsely rejecting the null hypothesis.




### Homogeneity of covariances assumption


box_m(nanopaper_stem[, "TI_star", drop = FALSE], nanopaper_stem$variety)






############### CLD ranking plot for lmer model ############### 


model.em <- emmeans(model,specs=~variety:section:HPH)

cld <- cld(model.em, alpha=0.05)

#print(cld$.group)

plot(cld) + 
  theme_bw() + 
  labs(x = "Estimated marginal mean", y = "Biomass samples") 
  

cld$ID <- str_c(cld$variety,' ',cld$section,' ',cld$HPH,' ')
cld$ID <- as.factor(cld$ID)



cld %>%
  mutate(ID = fct_reorder(ID, emmean,.fun='mean')) %>%
  ggplot(aes(x=emmean, y=ID))+
  geom_errorbar(aes(xmin=lower.CL,
                    xmax=upper.CL),
                size=5,
                width=0,
                colour = 'blue',
                alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  geom_text(aes(x=114,label=cld$.group),
            fontface=2,
            hjust = -0.5) +
  xlab("Mean Tensile Index (Nm/g)") +
  ylab("Sample ID") +
  theme_bw()









# GGplot with all sections, varieties & energy levels

ggplot(nanopaper_samples) +
  aes(x = HPH.fac, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter() +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("HPH energy level") +
  ylab("Tensile Index (Nm/g)")


# Sorghum Landscape plot (all data points)

#Weight v Thickness
ggplot(nanopaper_samples) +
  aes(x = weight, y = thickness, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Weight (mg)") +
  ylab("Thickness (mm)")



#TI v Strain
ggplot(nanopaper_samples) +
  aes(x = strain_star, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Strain at break (%)") +
  ylab("Tensile Index (Nm/g)")


#TI v Toughness
ggplot(nanopaper_samples) +
  aes(x = toughness_star, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Toughness (MJ/m3))") +
  ylab("Tensile Index (Nm/g)")


#TI v Modulus
ggplot(nanopaper_samples) +
  aes(x = E_star, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Modulus (GPA)") +
  ylab("Tensile Index (Nm/g)")


#Density v E
ggplot(nanopaper_samples) +
  aes(x = density, y = E_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Density (g/cm3)") +
  ylab("Modulus (GPa)")


#Density v Tensile 
ggplot(nanopaper_samples) +
  aes(x = density, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Density (g/cm3)") +
  ylab("Tensile Index (Nm/g)")


#Density v Toughness
ggplot(nanopaper_samples) +
  aes(x = density, y = toughness_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Density (g/cm3)") +
  ylab("Toughness (MJ/m3)")


#Grammage v TI
ggplot(nanopaper_samples) +
  aes(x = grammage, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Grammage ()") +
  ylab("Tensile Index (Nm/g)")

#Grammage v Toughness
ggplot(nanopaper_samples) +
  aes(x = grammage, y = toughness_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Grammage ()") +
  ylab("Toughness (MJ/m3)")


#Grammage v E
ggplot(nanopaper_samples) +
  aes(x = grammage, y = E_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Grammage ()") +
  ylab("Young's Modulus (GPa)")

#Weight v TI
ggplot(nanopaper_samples) +
  aes(x = weight, y = E_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter(alpha = 0.4) +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Grammage ()") +
  ylab("Young's Modulus (GPa)")



# Sorghum Landscape plot (average points)

grouped = nanopaper_samples %>%
  group_by(variety, section, HPH) %>%
  mutate(toughness_star = first(na.omit(toughness_star))) %>%
  mutate(TI_star = first(na.omit(TI_star))) %>%
  mutate(strain_star = first(na.omit(strain_star))) %>%
  mutate(E_star = first(na.omit(E_star)))


grouped_na_TI <- grouped %>%
  group_by(variety, section, HPH) %>%
  summarise(TI_star = mean(TI_star))

grouped_na_strain <- grouped %>%
  group_by(variety, section, HPH) %>%
  summarise(strain_star = mean(strain_star))

grouped_na_toughness <- grouped %>%
  group_by(variety, section, HPH) %>%
  summarise(toughness_star = mean(toughness_star))

grouped_na_E <- grouped %>%
  group_by(variety, section, HPH) %>%
  summarise(E_star = mean(E_star))

grouped_na_TI$strain_star <- grouped_na_strain$strain_star
grouped_na_TI$toughness_star <- grouped_na_toughness$toughness_star
grouped_na_TI$E_star <- grouped_na_E$E_star

grouped_na_TI$HPH <- factor(grouped_na_TI$HPH, levels = c("L", "M", "H"))

# Average TI vs strain  
ggplot(grouped_na_TI) +
  aes(x = strain_star, y = TI_star, color = variety, shape = section, size = HPH) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_jitter() +
  theme(legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Strain at break (%)") +
  ylab("Tensile Index (Nm/g)")



# Combined Sorghum Landscape plots

grouped_na_overlay <- grouped_na_TI

colnames(grouped_na_overlay)[1] <- "variety_ave"
colnames(grouped_na_overlay)[2] <- "section_ave"
colnames(grouped_na_overlay)[3] <- "HPH_ave"
colnames(grouped_na_overlay)[4] <- "TI_ave"
colnames(grouped_na_overlay)[5] <- "strain_ave"
colnames(grouped_na_overlay)[6] <- "toughness_ave"
colnames(grouped_na_overlay)[7] <- "E_ave"


grouped_na_overlay$ID <- str_c(grouped_na_overlay$variety_ave,'-',grouped_na_overlay$section_ave,'-',grouped_na_overlay$HPH_ave)
nanopaper_samples$ID <- str_c(nanopaper_samples$variety,'-',nanopaper_samples$section,'-',nanopaper_samples$HPH)

overlay <- dplyr::left_join(nanopaper_samples,grouped_na_overlay, by="ID")

HPH <- "HPH Energy Level"

#TI v Strain (overlay)
ggplot(overlay) +
  aes(x = strain_star, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  labs(color = "Variety", shape = "Section", size = "HPH Energy Level") +
  theme(axis.title=element_text(size=14,face="bold"), legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Strain at break (%)") +
  ylab("Tensile Index (Nm/g)") + 
  geom_jitter(alpha = 0.09) + 
  geom_point(mapping = aes(x = strain_ave, y = TI_ave))


#TI v Toughness (overlay)
ggplot(overlay) +
  aes(x = toughness_star, y = TI_star, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  labs(color = "Variety", shape = "Section", size = "HPH Energy Level") +
  theme(axis.title=element_text(size=14,face="bold"), legend.position = "right", panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Toughness (MJ/m3)") +
  ylab("Tensile Index (Nm/g)") + 
  geom_jitter(alpha = 0.09) + 
  geom_point(mapping = aes(x = toughness_ave, y = TI_ave))  










