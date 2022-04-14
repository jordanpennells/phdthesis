## Figures for Sorghum2020 publication

library(ggplot2)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(rstatix)
library(broom)
library(effects)
library(lme4)
library(Rcpp)
library(dplyr)
library(tidyverse)
library(multcomp)
library(report)
library(agricolae)


## Figure 1 - Sorghum Section Breakdown


##  mass_breakdown <- read_csv("mass breakdown.csv")

##  #mass_breakdown %>%
##    #bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))



##  bdtest <- read_csv("bdtest.csv")

##  ggplot(bdtest, aes(x=plant, y=leaf, fill=tiller))


## Analysis of variance

HPH.fac <- factor(nanopaper$HPH, levels = c("L", "M", "H"))


TI.aov <- aov(TI_star ~ variety + section + HPH.fac, data=nanopaper)
anova(TI.aov)
plot(allEffects(TI.aov))

strain.aov <- aov(strain_star ~ variety + section + HPH.fac, data=nanopaper)
anova(strain.aov)
plot(allEffects(strain.aov))

Ut.aov <- aov(toughness_star ~ variety + section + HPH.fac, data=nanopaper)
anova(Ut.aov)
plot(allEffects(Ut.aov))

E.aov <- aov(modulus_star ~ variety + section + HPH.fac, data=nanopaper)
anova(E.aov)
plot(allEffects(E.aov))


cv.model(TI.aov)
cv.model(strain.aov)
cv.model(Ut.aov)
cv.model(E.aov)

TI.aov.m <- aov(TI_star ~ variety * section * HPH.fac, data=nanopaper)
anova(TI.aov.m)
plot(allEffects(TI.aov.m))

strain.aov.m <- aov(strain_star ~ variety * section * HPH.fac, data=nanopaper)
anova(strain.aov.m)
plot(allEffects(strain.aov.m))

Ut.aov.m <- aov(toughness_star ~ variety * section * HPH.fac, data=nanopaper)
anova(Ut.aov.m)
plot(allEffects(Ut.aov.m))

E.aov.m <- aov(modulus_star ~ variety * section * HPH.fac, data=nanopaper)
anova(E.aov.m)
plot(allEffects(E.aov.m))


#Tukey test - Post Hoc ANOVA Test

TukeyHSD(TI.aov)
plot(TukeyHSD(TI.aov))


#Hsu MCB test - Post Hoc ANOVA Test

lmod <- lm(TI_star ~ variety+section+HPH, data = nanopaper)
aov(lmod)
report(lmod) 


nanopaper_samples$strip <- as.factor(nanopaper_samples$strip)
nanopaper_samples$nanopaper <- as.factor(nanopaper_samples$nanopaper)

model <- lmer(TI_star ~ HPH + variety + section + (section|variety) + (1|nanopaper), data = nanopaper_samples)

# model <- lmer(Strength ~ energy + variety + section + (section|variety) + (1|nanopaper), data = nanopaper_samples)


ggplot(nanopaper_samples,aes(x=section,y=TI_star,col=variety)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~variety) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=variety,y=TI_star,col=variety)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~nanopaper) +
  theme_bw()

ggplot(nanopaper_samples,aes(x=variety,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=section,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=nanopaper,y=TI_star,col=nanopaper)) + 
  geom_jitter() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=strip,y=TI_star,col=variety)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~variety) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=strip,y=TI_star,col=section)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~section) +
  theme_bw()


ggplot(nanopaper_samples,aes(x=strip,y=TI_star,col=HPH)) + 
  geom_point() + 
  geom_boxplot(alpha=0.2) + 
  facet_wrap(~HPH) +
  theme_bw()




## 3 Way ANOVA w/ replication


np_TI_summary <- nanopaper %>%
  group_by(variety, section, HPH) %>%
  get_summary_stats(TI_star, type = "mean_sd")

np_TI_summary_full <- nanopaper %>%
  group_by(variety, section, HPH) %>%
  get_summary_stats(TI_star, type = "full")


nanopaper_TI <- ggboxplot(nanopaper,x="HPH",y="TI_star",color="section",facet.by="variety",add="jitter", xlab="HPH Energy Level", ylab="Tensile Index (Nm/g)")
nanopaper_TI

nanopaper_UT <- ggboxplot(nanopaper,x="HPH",y="toughness_star",color="section",facet.by="variety",add="jitter", xlab="HPH Energy Level", ylab="Toughness (MJ/m3)")
nanopaper_UT

nanopaper_E <- ggboxplot(nanopaper,x="HPH",y="modulus_star",color="section",facet.by="variety",add="jitter", xlab="HPH Energy Level", ylab="Young's Modulus (GPa)")
nanopaper_E

      #Check normality assumption
qq  <- lm(TI_star ~ variety*section*HPH, data = nanopaper)
ggqqplot(residuals(qq))
shapiro.test(residuals(qq))
      #if the result is non-significant, assume normality

      #Check normality assumption by groups
nanopaper_normtest <- nanopaper %>%
  group_by(variety, section, HPH) %>%
  shapiro_test(TI_star)
      #all factors are normally distributed (p>0.03) except for Sugargraze-Stem(<1m)-L

#ggqqplot(nanopaper, "TI_star", ggtheme = theme_bw()) +
 # facet_grid(variety + section ~ HPH, labeller = "label_both")

      #Homogeneity of variance assumption
nanopaper %>% levene_test(TI_star ~ variety*section*HPH)
      #The Levene’s test is significant (p < 0.05), therefore we cannot assume the homogeneity of variances in the different groups.
      #The assumption of homogeneity of variance means that the level of variance for a particular variable is constant across the sample.
      # In ANOVA, when homogeneity of variance is violated there is a greater probability of falsely rejecting the null hypothesis.


threeway.aov_TI <- nanopaper %>% anova_test(TI_star ~ variety*section*HPH)
threeway.aov_TI
      #There was a statistically significant three-way interaction between variety, section, and HPH level
      #Continue with a simple two-way interaction ANOVA between variety and section


#Re-do Threeway ANOVA w/ stem sections combined

nanopaper_ANOVA <- read_csv("nanopaper_ANOVA.csv")

nanopaper_samples <- read_csv("nanopaper_samples.csv")

#nanopaper_ANOVA <- nanopaper_ANOVA %>% mutate_all(funs(str_replace(., "<1m", "")))
#nanopaper_ANOVA <- nanopaper_ANOVA %>% mutate_all(funs(str_replace(., ">1m", "")))


#nanopaper_ANOVA %>% anova_test(TI_star ~ variety + HPH + variety/section + variety:HPH + variety/section:HPH)

summary(aov(TI_star ~ variety + HPH + variety/section + variety:HPH + variety/section:HPH, nanopaper_ANOVA))
















##Simple two-way interaction ANOVA
#twoway_model <- lm(TI_star~variety*section*HPH, data=nanopaper)
#nanopaper %>%
  #group_by(HPH) %>%
  #anova_test(TI_star ~ variety*section, error = twoway_model)
      #For M/H energy, this result suggests that the effect of section on TI depends on the variety, but doesn't depend on variety for L energy.








## Boxplots

nanopaper$HPH <- factor(nanopaper$HPH, levels = c("L", "M", "H"))

bxp_TI_section <- ggboxplot(nanopaper, x = "section", y = "TI_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_TI_section

bxp_TI_variety <- ggboxplot(nanopaper, x = "variety", y = "TI_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_TI_variety

bxp_TI_variety2 <- ggplot(nanopaper, aes(x=variety, y=TI_star, colour = variety, fill=HPH)) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000")) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  geom_boxplot(alpha=0.3) +
  geom_point(position=position_jitterdodge(), aes(alpha=HPH)) +
  theme(axis.title=element_text(size=14,face="bold"), legend.position = "right", panel.background = element_rect(fill = "white",color="black")) +
  labs(color = "Variety", fill = "HPH Energy Level") +
  xlab("Sorghum Variety") +
  ylab("Tensile Index (Nm/g)") +
  guides(alpha="none")
bxp_TI_variety2

bxp_UT_section <- ggboxplot(nanopaper, x = "section", y = "toughness_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_UT_section

bxp_UT_variety <- ggboxplot(nanopaper, x = "variety", y = "toughness_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_UT_variety

bxp_UT_variety2 <- ggplot(nanopaper, aes(x=variety, y=toughness_star, colour = variety, fill=HPH)) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000")) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  geom_boxplot(alpha=0.3) +
  geom_point(position=position_jitterdodge(), aes(alpha=HPH)) +
  theme(axis.title=element_text(size=14,face="bold"), legend.position = "right", panel.background = element_rect(fill = "white",color="black")) +
  labs(color = "Variety", fill = "HPH Energy Level") +
  xlab("Sorghum Variety") +
  ylab("Toughness (MJ/m3)") +
  guides(alpha="none")
bxp_UT_variety2

bxp_E_section <- ggboxplot(nanopaper, x = "section", y = "modulus_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_E_section

bxp_E_variety <- ggboxplot(nanopaper, x = "variety", y = "modulus_star",color = "HPH", palette = c("#CCCCCC","#666666","#000000"),add="jitter")
bxp_E_variety

bxp_E_variety2 <- ggplot(nanopaper, aes(x=variety, y=modulus_star, colour = variety, fill=HPH)) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000")) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  geom_boxplot(alpha=0.3) +
  geom_point(position=position_jitterdodge(), aes(alpha=HPH)) +
  theme(axis.title=element_text(size=14,face="bold"), legend.position = "right", panel.background = element_rect(fill = "white",color="black")) +
  labs(color = "Variety", fill = "HPH Energy Level") +
  xlab("Sorghum Variety") +
  ylab("Young's Modulus (GPa)") +
  guides(alpha="none")
bxp_E_variety2


bxp_E_variety_violin <- ggplot(nanopaper, aes(x = variety.fac, y = modulus_star, fill = HPH)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000"))
bxp_E_variety_violin


bxp_TI_variety_violin <- ggplot(nanopaper, aes(x = HPH, y = TI_star, fill = HPH)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000")) +
  facet_wrap(~variety)
bxp_TI_variety_violin


bxp_TI_section_violin <- ggplot(nanopaper, aes(x = HPH, y = TI_star, fill = HPH)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=c("#CCCCCC","#666666","#000000")) +
  facet_wrap(~section)
bxp_TI_section_violin




##Linear Model

lm.model <- lm(TI_star ~ variety + section + HPH, data=nanopaper)
summary(lm.model)






