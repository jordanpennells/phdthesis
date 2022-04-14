

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
library(fmsb)





morfi_biomass <- read_csv("morfi_biomass.csv")

morfi_cnf<- read_csv("morfi_cnf.csv")




################################################################################

###Fibre Evolution Visualisation - Fibre Width

morfi_cnf$ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section,'-',morfi_cnf$HPH)

morfi_cnf %>% group_by(ID) %>% dplyr::summarise(n=n(),
                                         se=sd(fibre_W,na.rm=TRUE)/sqrt(n),
                                         mean=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean'))


  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3) +
  geom_point(stat="identity", size=2.5) +
  xlab("Mean Fibre Width (um)") +
  ylab("Sample ID") +
  theme_bw()



###Fibre Evolution Visualisation - Fibre Width - CLD plot

v1 <- max(abs(PCA.test$rotation[,1]),na.rm=TRUE)
v2 <- max(abs(PCA.test$rotation[,2]),na.rm=TRUE)
v3 <- max(abs(PCA.test$rotation[,3]),na.rm=TRUE)
v4 <- max(abs(PCA.test$rotation[,4]),na.rm=TRUE)

var1 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,1])==v1]
var2 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,2])==v2]
var3 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,3])==v3]
var4 <- rownames(PCA.test$rotation)[abs(PCA.test$rotation[,4])==v4]

aov.var1 <- aov(as.formula(paste(eval(var1),"HPH*Variety/Section",sep="~")), data = morfi_cnf)
MM.aov.var1 <- emmeans(aov.var1,specs=~HPH:Variety:Section)
cld.var1 <- cld(MM.aov.var1,reversed=TRUE,Letters=letters)

cld.var1$ID <- str_c(cld.var1$Variety,'-',cld.var1$Section,'-',cld.var1$HPH)
cld.var1$ID <- as.factor(cld.var1$ID)

morfi_cnf_emmeans <- merge(morfi_cnf, cld.var1)

morfi_cnf_emmeans %>% group_by(ID) %>% dplyr::summarise(n=n(),
                                                se=sd(fibre_W,na.rm=TRUE)/sqrt(n),
                                                mean=mean(fibre_W,na.rm=TRUE),
                                                .group=.group) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID= fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, tidytext::reorder_within(ID, mean, .group))) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3) +
  geom_point(stat="identity", size=2.5) +
  xlab("Mean Fibre Width (um)") +
  ylab("Sample ID") +
  geom_text(aes(x=30,label=morfi_cnf_emmeans$.group),fontface=2,hjust = 1) +
  theme_bw()










#Fibre Evolution Visualisation - Fibre Width - % Difference

morfi_cnf$biomass_ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section)
morfi_biomass$biomass_ID <- str_c(morfi_biomass$Variety,'-',morfi_biomass$Section)

biomass_stat <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n_biomass=n(),se_biomass=sd(fibre_W,na.rm=TRUE)/sqrt(n_biomass),mean_biomass=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(lower_biomass=mean_biomass-qt(0.975,n_biomass-1)*se_biomass,upper_biomass=mean_biomass+qt(0.975,n_biomass-1)*se_biomass) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean_biomass,.fun='mean'))

morfi_stat <- morfi_cnf %>% group_by(biomass_ID,HPH) %>% dplyr::summarise(n=n(),se=sd(fibre_W,na.rm=TRUE)/sqrt(n),mean=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  pivot_wider(id_cols=biomass_ID, values_from = c(n,se,mean,lower,upper),names_from = HPH) %>%
  full_join(biomass_stat) 
morfi_stat <- morfi_stat[,c(1,19,18,20,21,9,6,12,15,10,7,13,16,8,5,11,14)]

diff <- morfi_stat %>% mutate(lower_diff_BL = lower_biomass - upper_L, higher_diff_BL = upper_biomass - lower_L, lower_diff_LM = lower_L - upper_M, higher_diff_LM = upper_L - lower_M, lower_diff_MH = lower_M - upper_H, higher_diff_MH = upper_M - lower_H, BL = (higher_diff_BL - lower_diff_BL)/2 + lower_diff_BL, LM = (higher_diff_LM - lower_diff_LM)/2 + lower_diff_LM, MH = (higher_diff_MH - lower_diff_MH)/2 + lower_diff_MH)

diff_BL <- diff %>% ggplot(aes(x=BL, fct_reorder(biomass_ID,BL))) +
  labs(title="Difference between Biomass to Low", y = "Biomass sample", x = "% difference") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_diff_BL,
                     xmax=higher_diff_BL),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-1,42) +  theme_bw() 

diff_LM <- diff %>% ggplot(aes(x=LM, fct_reorder(biomass_ID,LM))) +
  labs(title="Difference between Low to Medium", y = "Biomass sample", x = "% difference") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_diff_LM,
                     xmax=higher_diff_LM),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-1,42)+  theme_bw() 

diff_MH <- diff %>% ggplot(aes(x=MH, fct_reorder(biomass_ID,MH))) +
  labs(title="Difference between Medium to High", x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_diff_MH,
                     xmax=higher_diff_MH),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-1,42)+  theme_bw() 

grid.arrange(diff_BL,diff_LM,diff_MH)


#Fibre Evolution Visualisation - Fibre Width - BLMH Panels

width_B <- morfi_stat %>% ggplot(aes(x=mean_biomass, fct_reorder(biomass_ID,mean_biomass))) +
  labs(title="Biomass powder fibre width", y = "Biomass sample", x = "Mean fibre width_Biomass powder") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_biomass,
                     xmax=upper_biomass),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(10,55) + theme_bw()

width_L <- morfi_stat %>% ggplot(aes(x=mean_L, fct_reorder(biomass_ID,mean_L))) +
  labs(title="Low energy fibre width", y = "Biomass sample", x = "Mean fibre width_Low energy") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_L,
                     xmax=upper_L),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(10,55) + theme_bw()

width_M <- morfi_stat %>% ggplot(aes(x=mean_M, fct_reorder(biomass_ID,mean_M))) +
  labs(title="Medium energy fibre width", y = "Biomass sample", x = "Mean fibre width_Medium energy") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_M,
                     xmax=upper_M),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(10,55) + theme_bw()  

width_H <- morfi_stat %>% ggplot(aes(x=mean_H, fct_reorder(biomass_ID,mean_H))) +
  labs(title="High energy fibre width", y = "Biomass sample", x = "Mean fibre width_High energy") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_H,
                     xmax=upper_H),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(10,55) +  theme_bw()   


grid.arrange(width_B,width_L,width_M,width_H, ncol=1) 


################ Fibre Evolution Visualisation - Fibre Width - Reaction norm


##morfi_stat_int <- morfi_stat %>% pivot_longer(cols=2:17,values_to=c("mean","se","lower","upper"))

morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,HPH) %>% dplyr::summarise(n=n(),se=sd(fibre_W,na.rm=TRUE)/sqrt(n),mean=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fibre_W,na.rm=TRUE)/sqrt(n),mean=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of fibre width", x = "Processing Stage", y = "Fibre width (um)") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 




################################################################################


#Fibre Evolution Visualisation - Fine Content Length-weighted Length - LMH

morfi_cnf$ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section,'-',morfi_cnf$HPH)

fine_cont.L.L <- morfi_cnf %>% group_by(ID) %>% dplyr::summarise(n=n(),se=sd(fine_cont.L.L,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont.L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  xlab("Mean Fine Content (% in length-weighted length)") +
  ylab("Sample ID")
fine_cont.L.L






#Fibre Evolution Visualisation - Fine Content Length-weighted Length - % Difference

morfi_cnf$biomass_ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section)
colnames(morfi_cnf)[5] <- "type"
morfi_biomass$biomass_ID <- str_c(morfi_biomass$Variety,'-',morfi_biomass$Section)

biomass_stat <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n_biomass=n(),se_biomass=sd(fine_cont.L.L,na.rm=TRUE)/sqrt(n_biomass),mean_biomass=mean(fine_cont.L.L,na.rm=TRUE)) %>%
  mutate(lower_biomass=mean_biomass-qt(0.975,n_biomass-1)*se_biomass,upper_biomass=mean_biomass+qt(0.975,n_biomass-1)*se_biomass) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean_biomass,.fun='mean'))

morfi_stat <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fine_cont.L.L,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont.L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  pivot_wider(id_cols=biomass_ID, values_from = c(n,se,mean,lower,upper),names_from = type) %>%
  full_join(biomass_stat) 
morfi_stat <- morfi_stat[,c(1,19,18,20,21,9,6,12,15,10,7,13,16,8,5,11,14)]

diff <- morfi_stat %>% mutate(lower_diff_BL = lower_biomass - upper_L, higher_diff_BL = upper_biomass - lower_L, lower_diff_LM = lower_L - upper_M, higher_diff_LM = upper_L - lower_M, lower_diff_MH = lower_M - upper_H, higher_diff_MH = upper_M - lower_H, BL = (higher_diff_BL - lower_diff_BL)/2 + lower_diff_BL, LM = (higher_diff_LM - lower_diff_LM)/2 + lower_diff_LM, MH = (higher_diff_MH - lower_diff_MH)/2 + lower_diff_MH)

diff_BL <- diff %>% ggplot(aes(x=-BL, fct_reorder(biomass_ID,BL))) +
  labs(title="Difference between Biomass to Low",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_BL,
                     xmax=-higher_diff_BL),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-300,200)+  theme_bw() 

diff_LM <- diff %>% ggplot(aes(x=-LM, fct_reorder(biomass_ID,LM))) +
  labs(title="Difference between Low to Medium",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_LM,
                     xmax=-higher_diff_LM),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-100,100)+  theme_bw() 

diff_MH <- diff %>% ggplot(aes(x=-MH, fct_reorder(biomass_ID,MH))) +
  labs(title="Difference between Medium to High", x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_MH,
                     xmax=-higher_diff_MH),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-100,100)+  theme_bw() 

grid.arrange(diff_BL,diff_LM,diff_MH)



#Fibre Evolution Visualisation - Fine Content Length-weighted Length - BLMH Panels

length.length_B <- morfi_stat %>% ggplot(aes(x=mean_biomass, fct_reorder(biomass_ID,mean_biomass))) +
  labs(title="Biomass Fine Content (% in length-weighted length)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_biomass,
                     xmax=upper_biomass),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400)+  theme_bw() 

length.length_L <- morfi_stat %>% ggplot(aes(x=mean_L, fct_reorder(biomass_ID,mean_L))) +
  labs(title="Low Energy Fine Content (% in length-weighted length)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_L,
                     xmax=upper_L),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400)+  theme_bw() 

length.length_M <- morfi_stat %>% ggplot(aes(x=mean_M, fct_reorder(biomass_ID,mean_M))) +
  labs(title="Medium Energy Fine Content (% in length-weighted length)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_M,
                     xmax=upper_M),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400) +  theme_bw()  

length.length_H <- morfi_stat %>% ggplot(aes(x=mean_H, fct_reorder(biomass_ID,mean_H))) +
  labs(title="High Energy Fine Content (% in length-weighted length)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_H,
                     xmax=upper_H),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400) +  theme_bw()   


grid.arrange(length.length_B,length.length_L,length.length_M,length.length_H, ncol=1) 



################ Fibre Evolution Visualisation - Fine Content Length-weighted Length - Reaction norm



morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fine_cont.L.L,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont.L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fine_cont.L.L,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont.L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of Fine Content (% in length-weighted length)", x = "Processing Stage", y = "Fine Content (% in length-weighted length)") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 







################################################################################


#Fibre Evolution Visualisation - Fibre Length-weighted Length - LMH

morfi_cnf$ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section,'-',morfi_cnf$HPH)

fibre_L.L <- morfi_cnf %>% group_by(ID) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  xlab("Mean Fibre Length-weighted Length (um)") +
  ylab("Sample ID")
fibre_L.L


#Fibre Evolution Visualisation - Fibre Length-weighted Length - % Difference

morfi_cnf$biomass_ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section)
colnames(morfi_cnf)[5] <- "type"
morfi_biomass$biomass_ID <- str_c(morfi_biomass$Variety,'-',morfi_biomass$Section)

biomass_stat <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n_biomass=n(),se_biomass=sd(fibre_L.L,na.rm=TRUE)/sqrt(n_biomass),mean_biomass=mean(fibre_L.L,na.rm=TRUE)) %>%
  mutate(lower_biomass=mean_biomass-qt(0.975,n_biomass-1)*se_biomass,upper_biomass=mean_biomass+qt(0.975,n_biomass-1)*se_biomass) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean_biomass,.fun='mean'))

morfi_stat <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  pivot_wider(id_cols=biomass_ID, values_from = c(n,se,mean,lower,upper),names_from = type) %>%
  full_join(biomass_stat) 
morfi_stat <- morfi_stat[,c(1,19,18,20,21,9,6,12,15,10,7,13,16,8,5,11,14)]

diff <- morfi_stat %>% mutate(lower_diff_BL = lower_biomass - upper_L, higher_diff_BL = upper_biomass - lower_L, lower_diff_LM = lower_L - upper_M, higher_diff_LM = upper_L - lower_M, lower_diff_MH = lower_M - upper_H, higher_diff_MH = upper_M - lower_H, BL = (higher_diff_BL - lower_diff_BL)/2 + lower_diff_BL, LM = (higher_diff_LM - lower_diff_LM)/2 + lower_diff_LM, MH = (higher_diff_MH - lower_diff_MH)/2 + lower_diff_MH)

diff_BL <- diff %>% ggplot(aes(x=-BL, fct_reorder(biomass_ID,BL))) +
  labs(title="Difference between Biomass to Low",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_BL,
                     xmax=-higher_diff_BL),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-300,200)+  theme_bw() 

diff_LM <- diff %>% ggplot(aes(x=-LM, fct_reorder(biomass_ID,LM))) +
  labs(title="Difference between Low to Medium",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_LM,
                     xmax=-higher_diff_LM),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-100,100)+  theme_bw() 

diff_MH <- diff %>% ggplot(aes(x=-MH, fct_reorder(biomass_ID,MH))) +
  labs(title="Difference between Medium to High", x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_MH,
                     xmax=-higher_diff_MH),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-100,100)+  theme_bw() 

grid.arrange(diff_BL,diff_LM,diff_MH)


#Fibre Evolution Visualisation - Fibre Length-weighted Length - BLMH Panels

length.length_B <- morfi_stat %>% ggplot(aes(x=mean_biomass, fct_reorder(biomass_ID,mean_biomass))) +
  labs(title="Biomass fibre length-weighted length", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_biomass,
                     xmax=upper_biomass),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400)+  theme_bw() 

length.length_L <- morfi_stat %>% ggplot(aes(x=mean_L, fct_reorder(biomass_ID,mean_L))) +
  labs(title="Low fibre length-weighted length", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_L,
                     xmax=upper_L),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400)+  theme_bw() 

length.length_M <- morfi_stat %>% ggplot(aes(x=mean_M, fct_reorder(biomass_ID,mean_M))) +
  labs(title="Medium fibre length-weighted length", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_M,
                     xmax=upper_M),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400) +  theme_bw()  

length.length_H <- morfi_stat %>% ggplot(aes(x=mean_H, fct_reorder(biomass_ID,mean_H))) +
  labs(title="High fibre length-weighted length", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_H,
                     xmax=upper_H),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,400) +  theme_bw()   


grid.arrange(length.length_B,length.length_L,length.length_M,length.length_H, ncol=1) 











################ Fibre Evolution Visualisation - Fibre Length - Reaction norm



morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fibre_L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of fibre length", x = "Processing Stage", y = "Fibre length (um)") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 





################ Fibre Evolution Visualisation - Fibre Length weighted length - Reaction norm



morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of fibre length-weighted length", x = "Processing Stage", y = "Fibre length-weighted length (um)") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 











################################################################################


#Fibre Evolution Visualisation - Fibre Kink Number - LMH

morfi_cnf$ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section,'-',morfi_cnf$type)

fibre_kink.num <- morfi_cnf %>% group_by(ID) %>% dplyr::summarise(n=n(),se=sd(fibre_kink.num,na.rm=TRUE)/sqrt(n),mean=mean(fibre_kink.num,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  xlab("Mean Fibre Kink Number") +
  ylab("Sample ID")
fibre_kink.num


#Fibre Evolution Visualisation - Fibre Fibre Kink Number - % Difference

morfi_cnf$biomass_ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section)
colnames(morfi_cnf)[5] <- "type"
morfi_biomass$biomass_ID <- str_c(morfi_biomass$Variety,'-',morfi_biomass$Section)

biomass_stat <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n_biomass=n(),se_biomass=sd(fibre_kink.num,na.rm=TRUE)/sqrt(n_biomass),mean_biomass=mean(fibre_kink.num,na.rm=TRUE)) %>%
  mutate(lower_biomass=mean_biomass-qt(0.975,n_biomass-1)*se_biomass,upper_biomass=mean_biomass+qt(0.975,n_biomass-1)*se_biomass) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean_biomass,.fun='mean'))

morfi_stat <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_kink.num,na.rm=TRUE)/sqrt(n),mean=mean(fibre_kink.num,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  pivot_wider(id_cols=biomass_ID, values_from = c(n,se,mean,lower,upper),names_from = type) %>%
  full_join(biomass_stat) 
morfi_stat <- morfi_stat[,c(1,19,18,20,21,9,6,12,15,10,7,13,16,8,5,11,14)]

diff <- morfi_stat %>% mutate(lower_diff_BL = lower_biomass - upper_L, higher_diff_BL = upper_biomass - lower_L, lower_diff_LM = lower_L - upper_M, higher_diff_LM = upper_L - lower_M, lower_diff_MH = lower_M - upper_H, higher_diff_MH = upper_M - lower_H, BL = (higher_diff_BL - lower_diff_BL)/2 + lower_diff_BL, LM = (higher_diff_LM - lower_diff_LM)/2 + lower_diff_LM, MH = (higher_diff_MH - lower_diff_MH)/2 + lower_diff_MH)

diff_BL <- diff %>% ggplot(aes(x=-BL, fct_reorder(biomass_ID,BL))) +
  labs(title="Difference between Biomass to Low",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_BL,
                     xmax=-higher_diff_BL),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-.5,.5)+  theme_bw() 

diff_LM <- diff %>% ggplot(aes(x=-LM, fct_reorder(biomass_ID,LM))) +
  labs(title="Difference between Low to Medium",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_LM,
                     xmax=-higher_diff_LM),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-.5,.5)+  theme_bw() 

diff_MH <- diff %>% ggplot(aes(x=-MH, fct_reorder(biomass_ID,MH))) +
  labs(title="Difference between Medium to High", x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_MH,
                     xmax=-higher_diff_MH),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-.5,.5)+  theme_bw() 

grid.arrange(diff_BL,diff_LM,diff_MH)


#Fibre Evolution Visualisation - Fibre Kink Number - BLMH Panels

length.length_B <- morfi_stat %>% ggplot(aes(x=mean_biomass, fct_reorder(biomass_ID,mean_biomass))) +
  labs(title="Biomass Fibre Kink Number", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_biomass,
                     xmax=upper_biomass),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0.74,1.55)+  theme_bw() 

length.length_L <- morfi_stat %>% ggplot(aes(x=mean_L, fct_reorder(biomass_ID,mean_L))) +
  labs(title="Low Fibre Kink Number", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_L,
                     xmax=upper_L),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0.74,1.55)+  theme_bw() 

length.length_M <- morfi_stat %>% ggplot(aes(x=mean_M, fct_reorder(biomass_ID,mean_M))) +
  labs(title="Medium Fibre Kink Number", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_M,
                     xmax=upper_M),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0.74,1.55) +  theme_bw()  

length.length_H <- morfi_stat %>% ggplot(aes(x=mean_H, fct_reorder(biomass_ID,mean_H))) +
  labs(title="High Fibre Kink Number", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_H,
                     xmax=upper_H),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0.74,1.55) +  theme_bw()   


grid.arrange(length.length_B,length.length_L,length.length_M,length.length_H, ncol=1) 




################ Fibre Evolution Visualisation - Fibre Kink Number - Reaction norm



morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_kink.num,na.rm=TRUE)/sqrt(n),mean=mean(fibre_kink.num,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fibre_kink.num,na.rm=TRUE)/sqrt(n),mean=mean(fibre_kink.num,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of Fibre Kink Number", x = "Processing Stage", y = "Fibre Kink Number") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 




# 
# ################ Fibre Evolution Visualisation - Fibre Length weighted length - Reaction norm
# 
# 
# 
# morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
#   mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
# colnames(morfi_stat_int)[2] <- "type"
# 
# 
# biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fibre_L.L,na.rm=TRUE)/sqrt(n),mean=mean(fibre_L.L,na.rm=TRUE)) %>%
#   mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
#   mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
#   mutate(type = "Biomass")
# 
# 
# overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
# overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])
# 
# ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
#   geom_line() +
#   geom_point(size = 3) + 
#   labs(title="Evolution of fibre length-weighted length", x = "Processing Stage", y = "Fibre length-weighted length (um)") +
#   scale_shape_manual(values=c(15,17,16,1,10)) +
#   scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
#   theme_bw() 



################################################################################


#Fibre Evolution Visualisation - Fine Content  - LMH

morfi_cnf$ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section,'-',morfi_cnf$HPH)

fine_cont <- morfi_cnf %>% group_by(ID) %>% dplyr::summarise(n=n(),se=sd(fine_cont,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3)+
  geom_point(stat="identity", size=2.5)+
  xlab("Mean Fine Content (%)") +
  ylab("Sample ID")
fine_cont



#Fibre Evolution Visualisation - Fine Content - % Difference

morfi_cnf$biomass_ID <- str_c(morfi_cnf$Variety,'-',morfi_cnf$Section)
colnames(morfi_cnf)[5] <- "type"
morfi_biomass$biomass_ID <- str_c(morfi_biomass$Variety,'-',morfi_biomass$Section)

biomass_stat <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n_biomass=n(),se_biomass=sd(fine_cont,na.rm=TRUE)/sqrt(n_biomass),mean_biomass=mean(fine_cont,na.rm=TRUE)) %>%
  mutate(lower_biomass=mean_biomass-qt(0.975,n_biomass-1)*se_biomass,upper_biomass=mean_biomass+qt(0.975,n_biomass-1)*se_biomass) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean_biomass,.fun='mean'))

morfi_stat <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fine_cont,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  pivot_wider(id_cols=biomass_ID, values_from = c(n,se,mean,lower,upper),names_from = type) %>%
  full_join(biomass_stat) 
morfi_stat <- morfi_stat[,c(1,19,18,20,21,9,6,12,15,10,7,13,16,8,5,11,14)]

diff <- morfi_stat %>% mutate(lower_diff_BL = lower_biomass - upper_L, higher_diff_BL = upper_biomass - lower_L, lower_diff_LM = lower_L - upper_M, higher_diff_LM = upper_L - lower_M, lower_diff_MH = lower_M - upper_H, higher_diff_MH = upper_M - lower_H, BL = (higher_diff_BL - lower_diff_BL)/2 + lower_diff_BL, LM = (higher_diff_LM - lower_diff_LM)/2 + lower_diff_LM, MH = (higher_diff_MH - lower_diff_MH)/2 + lower_diff_MH)

diff_BL <- diff %>% ggplot(aes(x=-BL, fct_reorder(biomass_ID,BL))) +
  labs(title="Difference between Biomass to Low",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_BL,
                     xmax=-higher_diff_BL),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-400,600)+  theme_bw() 

diff_LM <- diff %>% ggplot(aes(x=-LM, fct_reorder(biomass_ID,LM))) +
  labs(title="Difference between Low to Medium",x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_LM,
                     xmax=-higher_diff_LM),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-400,600)+  theme_bw() 

diff_MH <- diff %>% ggplot(aes(x=-MH, fct_reorder(biomass_ID,MH))) +
  labs(title="Difference between Medium to High", x = "% difference", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=-lower_diff_MH,
                     xmax=-higher_diff_MH),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(-400,600)+  theme_bw() 

grid.arrange(diff_BL,diff_LM,diff_MH)



#Fibre Evolution Visualisation - Fine Content - BLMH Panels

length.length_B <- morfi_stat %>% ggplot(aes(x=mean_biomass, fct_reorder(biomass_ID,mean_biomass))) +
  labs(title="Biomass Fine Content (%)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_biomass,
                     xmax=upper_biomass),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,1100)+  theme_bw() 

length.length_L <- morfi_stat %>% ggplot(aes(x=mean_L, fct_reorder(biomass_ID,mean_L))) +
  labs(title="Low Energy Fine Content (%)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_L,
                     xmax=upper_L),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,1100)+  theme_bw() 

length.length_M <- morfi_stat %>% ggplot(aes(x=mean_M, fct_reorder(biomass_ID,mean_M))) +
  labs(title="Medium Energy Fine Content (%)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_M,
                     xmax=upper_M),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,1100) +  theme_bw()  

length.length_H <- morfi_stat %>% ggplot(aes(x=mean_H, fct_reorder(biomass_ID,mean_H))) +
  labs(title="High Energy Fine Content (%)", y = "Biomass sample") +
  geom_point(stat="identity", size=1.5) +
  geom_errorbarh(aes(xmin=lower_H,
                     xmax=upper_H),
                 colour = 'blue',
                 alpha=0.3)+
  xlim(0,1100) +  theme_bw()   


grid.arrange(length.length_B,length.length_L,length.length_M,length.length_H, ncol=1) 



################ Fibre Evolution Visualisation - Fine Content - Reaction norm



morfi_stat_int <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(n=n(),se=sd(fine_cont,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se)
colnames(morfi_stat_int)[2] <- "type"


biomass_stat_int <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(n=n(),se=sd(fine_cont,na.rm=TRUE)/sqrt(n),mean=mean(fine_cont,na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(biomass_ID = fct_reorder(biomass_ID, mean,.fun='mean')) %>%
  mutate(type = "Biomass")


overall_int <- bind_rows(biomass_stat_int,morfi_stat_int) %>% mutate(type=ordered(type,levels=c("Biomass","L","M","H")))
overall_int <- overall_int %>% mutate(Variety = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,1],section = str_split(overall_int$biomass_ID,"-",simplify = TRUE)[,2])

ggplot(overall_int,aes(x=type, y=mean, color=Variety, shape = section, group = biomass_ID)) +
  geom_line() +
  geom_point(size = 3) + 
  labs(title="Evolution of Fine Content (%)", x = "Processing Stage", y = "Fine Content (million/g pulp)") +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  theme_bw() 
































.#################### Fibre Evolution - Radar Plot ####################


### Fibre_L ###

morfi_evolution_fibre.L <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(mean_fibre_L=mean(fibre_L,na.rm=TRUE)) 


biomass_evolution_fibre.L <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(mean_fibre_L=mean(fibre_L,na.rm=TRUE)) %>%
  mutate(type = "Biomass")


### Fibre_W ###

morfi_evolution_fibre.W <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(mean_fibre_W=mean(fibre_W,na.rm=TRUE)) 


biomass_evolution_fibre.W <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(mean_fibre_W=mean(fibre_W,na.rm=TRUE)) %>%
  mutate(type = "Biomass")


### Fine_A ###

morfi_evolution_fine.A <- morfi_cnf %>% group_by(biomass_ID,type) %>% dplyr::summarise(mean_fine_n=mean(fine_n,na.rm=TRUE)) 


biomass_evolution_fine.A <- morfi_biomass %>% group_by(biomass_ID) %>% dplyr::summarise(mean_fine_n=mean(fine_n,na.rm=TRUE)) %>%
  mutate(type = "Biomass")





overall_evolution <- bind_rows(biomass_evolution_fibre.L,
                               biomass_evolution_fibre.W,
                               biomass_evolution_fine.A,
                               morfi_evolution_fibre.L,
                               morfi_evolution_fibre.W,
                               morfi_evolution_fine.A) %>% 
  mutate(type=ordered(type,levels=c("Biomass","L","M","H")))


overall_evolution <- overall_evolution %>% 
  mutate(Variety = str_split(overall_evolution$biomass_ID,"-",simplify = TRUE)[,1],
         section = str_split(overall_evolution$biomass_ID,"-",simplify = TRUE)[,2])








overall_evolution.test <- overall_evolution %>% dplyr::summarise(type,mean_fibre_W,mean_fibre_L)

radarchart(overall_evolution[c("mean_fibre_L","mean_fibre_W","mean_fine_n")])

radarchart(overall_evolution[c("Biomass","L","M","H")])










################ Fine:Fibre Ratio ################


morfi_cnf_FF <- read_csv("morfi_cnf_FineFibreRatio.csv")

morfi_cnf_FF$ID <- str_c(morfi_cnf_FF$Variety,'-',
                         morfi_cnf_FF$Section,'-',
                         morfi_cnf_FF$HPH)

morfi_cnf_FF$fine_fibre <- as.numeric(morfi_cnf_FF$fine_fibre)

morfi_cnf_FF$HPH.fac <- factor(morfi_cnf_FF$HPH, levels = c("L", "M", "H"))

morfi_cnf_FF %>% group_by(ID) %>% 
  dplyr::summarise(n=n(), 
                   se=sd(fine_fibre,na.rm=TRUE)/sqrt(n),
                   mean=mean(fine_fibre, na.rm=TRUE)) %>%
  mutate(lower=mean-qt(0.975,n-1)*se,upper=mean+qt(0.975,n-1)*se) %>%
  mutate(ID = fct_reorder(ID, mean,.fun='mean')) %>%
  ggplot(aes(x=mean, y=ID)) +
  geom_errorbarh(aes(xmin=lower,
                     xmax=upper),
                 colour = 'blue',
                 alpha=0.3) +
  geom_point(stat="identity", size=2.5) +
  xlab("Fine:Fibre Content Ratio") +
  ylab("Sample ID") +
  theme_bw()

ggplot(morfi_cnf_FF,aes(x=HPH.fac,y=fine_fibre,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~Section ~Stature, nrow=1) +
  theme_bw()




ggplot(morfi_cnf_FF,aes(x=HPH.fac,y=fine_cont,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~Section) +
  theme_bw()


ggplot(morfi_cnf_FF,aes(x=HPH.fac,y=fibre_cont,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~Section) +
  theme_bw()












ggplot(morfi_cnf_FF,aes(x=Section,y=fine_fibre,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~HPH.fac) +
  theme_bw() +
  ggtitle("Evolution of (Fine:Fibre) Ratio over Energy Series")


ggplot(morfi_cnf_FF,aes(x=Section,y=fine_cont,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~HPH.fac) +
  theme_bw() +
  ggtitle("Evolution of Fine Content over Energy Series")


ggplot(morfi_cnf_FF,aes(x=Section,y=fibre_cont,fill=Section)) + 
  geom_bar(stat="identity",position=position_dodge()) + 
  facet_wrap(~HPH.fac) +
  theme_bw() +
  ggtitle("Evolution of Fibre Content over Energy Series")
  



# morfi_cnf_FF_na <- morfi_cnf_FF[complete.cases(morfi_cnf_FF), ] 
# 
# morfi_cnf_FF2 <- subset(morfi_cnf_FF,select = -c(ID,Sample_name,Variety,Replicate,Section,HPH,Morfi.Rep))
# 
# morfi_cnf_FF2_PCA <- prcomp(~ ., data=morfi_cnf_FF2, na.action=na.omit, center = TRUE, scale=TRUE)
# 
# 
# summary(morfi_cnf_FF2_PCA)
# 
# fviz_eig(morfi_cnf_FF2_PCA)
# 
# ggbiplot(morfi_cnf_FF2_PCA)
# 
# morfi_cnf_FF2.section <- morfi_cnf_FF_na$Section
# 
# p_section <- ggbiplot(morfi_cnf_FF2_PCA,
#                       varname.size = 3,
#                       obs.scale = 0.5, 
#                       ellipse=TRUE, 
#                       ellipse.prob = 0.68, 
#                       alpha = 0.5, 
#                       groups = morfi_cnf_FF2.section) + 
#   scale_color_discrete(name = 'Sections') + 
#   theme(legend.direction = 'horizontal', 
#         legend.position = 'top', 
#         panel.background = element_blank())+
#   geom_hline(yintercept=0,colour = "grey80", size = 0.5, lty="dashed")+
#   geom_vline(xintercept=0,colour = "grey80", size = 0.5, lty="dashed")
# 
# p_section + theme(panel.grid.minor = element_blank(),
#                   panel.grid.major = element_blank(),
#                   panel.background = element_rect(fill="white",colour = "black", size = 2)) +
#   ggtitle("PCA: Nanopaper + Fibre-Water Interaction Metrics")











