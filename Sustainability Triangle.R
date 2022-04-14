
#Install libraries

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(ggmosaic)
library(nlme)
library(tidyr)
library(broom)



    #Preparing data

nanopaper_samples <- read_csv("nanopaper_samples.csv")

nanopaper_esr <- dplyr::select(nanopaper_samples,c(2,3,4,7:16))



    #Organising Tensile Index (TI) parameters

nanopaper_esr.TI <- dplyr::select(nanopaper_esr,c(1:3,13))
nanopaper_esr.TI <- na.omit(nanopaper_esr.TI)

nanopaper_esr.TI$HPH <- gsub("L",0,gsub("M",0.5,gsub("H",1,nanopaper_esr.TI$HPH)))

nanopaper_esr.TI$ID = str_c(nanopaper_esr.TI$variety,'-',nanopaper_esr.TI$section)

nanopaper_norm.TI <- nanopaper_esr.TI %>% group_by(ID,HPH) %>% 
  summarise(mean(TI_star))

nanopaper_esr.TI <- left_join(nanopaper_esr.TI,nanopaper_norm,by=c("ID","HPH"))

nanopaper_esr.TI$TI_norm <- with(nanopaper_esr.TI,(TI_star-min(`mean(TI_star)`))/(diff(range(`mean(TI_star)`))))


#nanopaper_esr.TI$TI_norm <- unlist(nanopaper_esr.TI$TI_norm)


    ##ESR Plotting - Tensile Index

ggplot(nanopaper_esr.TI) +
  aes(x = HPH, y = TI_norm, color = variety, shape = section) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point() + 
  stat_summary(aes(group=ID), fun=mean, geom="line", colour="green") +
  ylab("Tensile Index (normalised)")


    ##Sustainability Triangle - Tensile Index


fits <- lmList(TI_norm ~ as.numeric(HPH) | ID, data=nanopaper_esr.TI)

fits.coef <- coef(fits)

colnames(fits.coef) <- c("intercept","slope")

fits.coef$ID <- rownames(fits.coef)
IDsplit <- str_split(fits.coef$ID,"-",simplify = TRUE)
fits.coef$variety <- IDsplit[,1]
fits.coef$section <- IDsplit[,2]

fits.coef$angle <- atan(fits.coef$intercept/fits.coef$slope)
fits.coef$mag <- sqrt((fits.coef$intercept^2)+(fits.coef$slope^2))



#fits.coef$rotate.m <- cos(pi/4)*fits.coef$slope - sin(pi/4)*fits.coef$intercept
#fits.coef$rotate.c <- cos(pi/4)*fits.coef$intercept + sin(pi/4)*fits.coef$slope


#plot(fits.coef$rotate.m,fits.coef$rotate.c, xlim=c(-0.71,0.71),ylim=c(0,0.71),asp=1)

max(nanopaper_esr.TI$`mean(TI_star)`)

ggplot(fits.coef) +
  aes(x=slope,y=intercept, color = variety, shape = section) + 
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1, size = 1) +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size=12), axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.4), plot.subtitle = element_text(hjust = 0.4, size = 14), plot.caption = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", color = "Variety", shape = "Section", title = "Processing Sustainability Triangle", subtitle = "Tensile Index") +
  guides(color = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.5), shape = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.4))
  










#, caption = max(nanopaper_esr.TI$`mean(TI_star)`)  
################################################################################  

    #Organising Toughness (T) parameters

nanopaper_esr.T <- dplyr::select(nanopaper_esr,c(1:3,12))
nanopaper_esr.T <- na.omit(nanopaper_esr.T)

nanopaper_esr.T$HPH <- gsub("L",0,gsub("M",0.5,gsub("H",1,nanopaper_esr.T$HPH)))

nanopaper_esr.T$ID = str_c(nanopaper_esr.T$variety,'-',nanopaper_esr.T$section)

nanopaper_norm.T <- nanopaper_esr.T %>% group_by(ID,HPH) %>% 
  dplyr::summarise(mean(toughness_star))

nanopaper_esr.T <- left_join(nanopaper_esr.T,nanopaper_norm.T,by=c("ID","HPH"))

nanopaper_esr.T$T_norm <- with(nanopaper_esr.T,(toughness_star-min(`mean(toughness_star)`))/(diff(range(`mean(toughness_star)`))))

#nanopaper_esr.T$T_norm <- unlist(nanopaper_esr.T$T_norm)

    ##ESR Plotting - Toughness

ggplot(nanopaper_esr.T) +
  aes(x = HPH, y = T_norm, color = variety, shape = section) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point() + 
  ylab("Toughness (normalised)")


    ##Sustainability Triangle - Toughness

fits <- lmList(T_norm ~ as.numeric(HPH) | ID, data=nanopaper_esr.T)
fits.coef <- coef(fits)
colnames(fits.coef) <- c("intercept","slope")

fits.coef$ID <- rownames(fits.coef)
IDsplit <- str_split(fits.coef$ID,"-",simplify = TRUE)
fits.coef$variety <- IDsplit[,1]
fits.coef$section <- IDsplit[,2]


ggplot(fits.coef) +
  aes(x=slope,y=intercept, color = variety, shape = section) + 
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1, size  = 1) +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.4), plot.subtitle = element_text(hjust = 0.4, size = 14), plot.caption = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", color = "Variety", shape = "Section", title = "Processing Sustainability Triangle", subtitle = "Toughness") +
  guides(color = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.5), shape = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.4))


################################################################################  

    #Organising Modulus (E) parameters

nanopaper_esr.E <- dplyr::select(nanopaper_esr,c(1:3,11))
nanopaper_esr.E <- na.omit(nanopaper_esr.E)

nanopaper_esr.E$HPH <- gsub("L",0,gsub("M",0.5,gsub("H",1,nanopaper_esr.E$HPH)))

nanopaper_esr.E$ID = str_c(nanopaper_esr.E$variety,'-',nanopaper_esr.E$section)

nanopaper_norm.E <- nanopaper_esr.E %>% group_by(ID,HPH) %>% 
  dplyr::summarise(mean(E_star))

nanopaper_esr.E <- left_join(nanopaper_esr.E,nanopaper_norm.E,by=c("ID","HPH"))

nanopaper_esr.E$E_norm <- with(nanopaper_esr.E,(E_star-min(`mean(E_star)`))/(diff(range(`mean(E_star)`))))



    ##ESR Plotting - Modulus

ggplot(nanopaper_esr.E) +
  aes(x = HPH, y = E_norm, color = variety, shape = section) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point() + 
  ylab("Young's Modulus (normalised)")


    ##Sustainability Triangle - Modulus

fits <- lmList(E_norm ~ as.numeric(HPH) | ID, data=nanopaper_esr.E)
fits.coef <- coef(fits)
colnames(fits.coef) <- c("intercept","slope")

fits.coef$ID <- rownames(fits.coef)
IDsplit <- str_split(fits.coef$ID,"-",simplify = TRUE)
fits.coef$variety <- IDsplit[,1]
fits.coef$section <- IDsplit[,2]


ggplot(fits.coef) +
  aes(x=slope,y=intercept, color = variety, shape = section) + 
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1,  size=1) +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.4), plot.subtitle = element_text(hjust = 0.4, size = 14), plot.caption = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", color = "Variety", shape = "Section", title = "Processing Sustainability Triangle", subtitle = "Young's Modulus") +
  guides(color = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.5), shape = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.4))


















################################################################################  

#Organising Hydrogel Injectability parameters

hydrogel_esr <- read_csv("hydrogel_samples.csv")

hydrogel_esr$HPH <- gsub("L",0,gsub("M",0.5,gsub("H",1,hydrogel_esr$HPH)))

hydrogel_esr$ID = str_c(hydrogel_esr$variety,'-',hydrogel_esr$section)

hydrogel_norm <- hydrogel_esr %>% group_by(ID,HPH) %>% 
  dplyr::summarise(mean(glide))

hydrogel_esr <- left_join(hydrogel_esr,hydrogel_norm,by=c("ID","HPH"))

hydrogel_esr$glide_norm <- with(hydrogel_esr,(glide-min(`mean(glide)`))/(diff(range(`mean(glide)`))))




###### Organising Hydrogel Injectability parameters - Low/Medium

hydrogel_esr_LM <- read_csv("hydrogel_samples_LM.csv")

hydrogel_esr_LM$HPH <- gsub("L",0,gsub("M",1,hydrogel_esr_LM$HPH))

hydrogel_esr_LM$ID = str_c(hydrogel_esr_LM$variety,'-',hydrogel_esr_LM$section)

hydrogel_norm_LM <- hydrogel_esr_LM %>% group_by(ID,HPH) %>% 
  dplyr::summarise(mean(glide))

hydrogel_esr_LM <- left_join(hydrogel_esr_LM,hydrogel_norm_LM,by=c("ID","HPH"))

hydrogel_esr_LM$glide_norm <- with(hydrogel_esr_LM,(glide-min(`mean(glide)`))/(diff(range(`mean(glide)`))))


##ESR Plotting - Hydrogel

ggplot(hydrogel_esr_LM) +
  aes(x = HPH, y = glide_norm, color = variety, shape = section) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point() + 
  stat_summary(aes(group=ID), fun.y=mean, geom="line", colour="purple") +
  ylab("Injectability Glide Force (normalised)")


##Sustainability Triangle - Hydrogel

hydrogel_fits <- lmList(glide_norm ~ as.numeric(HPH) | ID, data=hydrogel_esr_LM)
hydrogel_fits.coef <- coef(hydrogel_fits)
colnames(hydrogel_fits.coef) <- c("intercept","slope")

hydrogel_fits.coef$ID <- rownames(hydrogel_fits.coef)
IDsplit <- str_split(hydrogel_fits.coef$ID,"-",simplify = TRUE)
hydrogel_fits.coef$variety <- IDsplit[,1]
hydrogel_fits.coef$section <- IDsplit[,2]


hydrogel_fits.coef$variety <- paste0("Hydrogel-", hydrogel_fits.coef$variety)


#fits.coef$slope[1] <- (hydrogel_esr[2,7] - hydrogel_esr[1,7]) / 0.5
#fits.coef$slope[2] <- (hydrogel_esr[5,7] - hydrogel_esr[4,7]) / 0.5
#fits.coef$slope[3] <- (hydrogel_esr[8,7] - hydrogel_esr[7,7]) / 0.5
#fits.coef$slope[4] <- (hydrogel_esr[11,7] - hydrogel_esr[10,7]) / 0.5
#fits.coef$slope <- unlist(fits.coef$slope)

ggplot(hydrogel_fits.coef) +
  aes(x=slope,y=intercept, color = variety, shape = section) + 
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1,  size=1) +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.4), plot.subtitle = element_text(hjust = 0.4, size = 14), plot.caption = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", color = "Variety", shape = "Section", title = "Processing Sustainability Triangle", subtitle = "Hydrogel Injectability Force") +
  guides(color = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.5), shape = guide_legend(title.position = "top", label.position = "bottom", title.vjust = 0.8, title.hjust = 0.4))





################################################################################  

### Combined Hydrogel + Nanopaper Triangle ###


combine.fits.coef <- fits.coef[,1:5]


nanopaper.hydrogel_fits.coef <- data.frame(rbind(as.matrix(combine.fits.coef), 
                                                 as.matrix(hydrogel_fits.coef)))


nanopaper.hydrogel_fits.coef$intercept <- as.numeric(nanopaper.hydrogel_fits.coef$intercept)
nanopaper.hydrogel_fits.coef$slope <- as.numeric(nanopaper.hydrogel_fits.coef$slope)

nanopaper.hydrogel_fits.coef$section <- as.factor(nanopaper.hydrogel_fits.coef$section)
nanopaper.hydrogel_fits.coef$variety <- factor(nanopaper.hydrogel_fits.coef$variety, 
                                               levels = c("Graingrass",
                                                          "GreenleafBMR",
                                                          "Sugargraze",
                                                          "Yemen",
                                                          "Hydrogel-Sugargraze"))

nanopaper.hydrogel_fits.coef$PointSize <- c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,10,10,10,10)



ggplot(nanopaper.hydrogel_fits.coef) +
  aes(x=slope,y=intercept, color = variety, shape = section) + 
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange","red")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(aes(size = PointSize)) +
  scale_size_binned_area(
    limits = c(-1, 10),
    breaks = c(5, 7)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 11) +
  geom_abline(slope = -1, intercept = 1, size = 1) +
  theme(aspect.ratio = 1, 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold", size = 14), 
        legend.text = element_text(size=12), 
        axis.title = element_text(face = "bold", size = 14), 
        plot.title = element_text(face = "bold", size = 20, hjust = 0.4), 
        plot.subtitle = element_text(hjust = 0.4, size = 14), 
        plot.caption = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1)) +
  labs(x = "Slope", y = "Intercept", 
       color = "Variety", shape = "Section", 
       title = "Processing Sustainability Triangle", 
       subtitle = "Nanopaper Tensile Index vs. Hydrogel Injectability") +
  guides(size = FALSE,color = guide_legend(title.position = "top", 
                              label.position = "bottom", 
                              title.vjust = 0.8, 
                              title.hjust = 0.5), 
         shape = guide_legend(title.position = "top", 
                              label.position = "bottom", 
                              title.vjust = 0.8, 
                              title.hjust = 0.4))




############################









