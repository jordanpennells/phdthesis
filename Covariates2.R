

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(olsrr)
library(scales)
library(csv)
theme_set(theme_classic())


################################# Biochemical Covariate ###########################



## Analysing covariate data - Biochemical composition

biochem_comp2 <- read_csv("Biochem comp2.csv")


## Generating new variable based on wt% basis and CHL basis

biochem_comp2$C_wt <- biochem_comp2$Glucan
biochem_comp2$H_wt <- biochem_comp2$Xylan + biochem_comp2$Mannan + biochem_comp2$Arabinan + biochem_comp2$Galactan + biochem_comp2$Rhamnan
biochem_comp2$L_wt <- biochem_comp2$Klason_Lignin + biochem_comp2$AS_Lignin
biochem_comp2$E_wt <- biochem_comp2$Extractives
biochem_comp2$A_wt <- biochem_comp2$Ash + biochem_comp2$Acid_Insoluble_Ash

biochem_comp2$C_CHL <- biochem_comp2$C_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt)
biochem_comp2$H_CHL <- biochem_comp2$H_wt/(biochem_comp2$H_wt + biochem_comp2$H_wt + biochem_comp2$L_wt)
biochem_comp2$L_CHL <- biochem_comp2$L_wt/(biochem_comp2$L_wt + biochem_comp2$H_wt + biochem_comp2$L_wt)




biochem_comp2$biochem.ID <- paste(biochem_comp2$Variety, biochem_comp2$Section, sep = "-")

fits.coef$biochem.ID <- rownames(fits.coef)

fits.coef$mag <- sqrt((fits.coef$intercept^2)+(fits.coef$slope^2))

biochem_fit.coef2 <- merge(biochem_comp2,fits.coef,by="biochem.ID")

biochem_fit.coef2$mc <- biochem_fit.coef2$slope + biochem_fit.coef2$intercept


### Biochemical composition - wt % basis

  #Finding best parameters to include in MAGNITUDE model

biochem_model2_mag <- lm(mag ~ C_wt + H_wt + L_wt + E_wt + A_wt + 0 , data = biochem_fit.coef2)

ols_step_all_possible(biochem_model2_mag)

ols_step_best_subset(biochem_model2_mag)

biochem_model2_mag <- lm(mag ~ C_wt + L_wt + E_wt + A_wt + 0, data = biochem_fit.coef2)
summary(biochem_model2_mag)


  #Finding best parameters to include in MC model

biochem_model2_mc <- lm(mc ~ C_wt + H_wt + L_wt + E_wt + A_wt + 0, data = biochem_fit.coef2)

ols_step_all_possible(biochem_model2_mc)

ols_step_best_subset(biochem_model2_mc)


biochem_model2_mc <- lm(mc ~ L_wt + E_wt + A_wt + 0, data = biochem_fit.coef2)
summary(biochem_model2_mc)


  #Finding best parameters to include in INTERCEPT model
  
biochem_model2_int <- lm(intercept ~ C_wt + H_wt + L_wt + E_wt + A_wt + 0, data = biochem_fit.coef2)

ols_step_all_possible(biochem_model2_int)

ols_step_best_subset(biochem_model2_int)


biochem_model2_int <- lm(intercept ~ A_wt + 0, data = biochem_fit.coef2)
summary(biochem_model2_int)



### Biochemical composition - CHL % basis

  #Finding best parameters to include in CHL-MAGNITUDE model

biochem_CHL_model_mag <- lm(mag ~ C_CHL + H_CHL + L_CHL + 0 , data = biochem_fit.coef2)

ols_step_all_possible(biochem_CHL_model_mag)

ols_step_best_subset(biochem_CHL_model_mag)

biochem_CHL_model2_mag <- lm(mag ~ C_CHL + L_CHL + E_CHL + A_CHL + 0, data = biochem_fit.coef2)
summary(biochem_CHL_model2_mag)



  #Finding best parameters to include in MC model

biochem_CHL_model_mc <- lm(mc ~ C_CHL + H_CHL + L_CHL + E_CHL + A_CHL + 0 , data = biochem_fit.coef2)

ols_step_all_possible(biochem_CHL_model_mc)

ols_step_best_subset(biochem_CHL_model_mc)

biochem_CHL_model2_mc <- lm(mag ~ L_CHL + E_CHL + A_CHL + 0, data = biochem_fit.coef2)
summary(biochem_CHL_model2_mc)

qqnorm(biochem_CHL_model2_mc$residuals)


### Biochemical composition - CHLEA ratios

biochem_comp2$C_CHLEA <- biochem_comp2$C_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt + biochem_comp2$E_wt + biochem_comp2$A_wt)
biochem_comp2$H_CHLEa <- biochem_comp2$H_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt + biochem_comp2$E_wt + biochem_comp2$A_wt)
biochem_comp2$L_CHLEA <- biochem_comp2$L_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt + biochem_comp2$E_wt + biochem_comp2$A_wt)
biochem_comp2$E_CHLEA <- biochem_comp2$E_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt + biochem_comp2$E_wt + biochem_comp2$A_wt)
biochem_comp2$A_CHLEA <- biochem_comp2$A_wt/(biochem_comp2$C_wt + biochem_comp2$H_wt + biochem_comp2$L_wt + biochem_comp2$E_wt + biochem_comp2$A_wt)

biochem_comp2$C.H <- biochem_comp2$C_CHLEA/biochem_comp2$H_CHLEA
biochem_comp2$C.L <- biochem_comp2$C_CHLEA/biochem_comp2$L_CHLEA
biochem_comp2$C.E <- biochem_comp2$C_CHLEA/biochem_comp2$E_CHLEA
biochem_comp2$C.A <- biochem_comp2$C_CHLEA/biochem_comp2$A_CHLEA
biochem_comp2$H.L <- biochem_comp2$H_CHLEA/biochem_comp2$L_CHLEA
biochem_comp2$H.E <- biochem_comp2$H_CHLEA/biochem_comp2$E_CHLEA
biochem_comp2$H.A <- biochem_comp2$H_CHLEA/biochem_comp2$A_CHLEA
biochem_comp2$L.E <- biochem_comp2$L_CHLEA/biochem_comp2$E_CHLEA
biochem_comp2$L.A <- biochem_comp2$L_CHLEA/biochem_comp2$A_CHLEA
biochem_comp2$E.A <- biochem_comp2$E_CHLEA/biochem_comp2$A_CHLEA

biochem_CHL.ratio_model_mc <- lm(mc ~ C.H + C.L + H.L + 0 , data = biochem_fit.coef2)

ols_step_all_possible(biochem_CHL_model_mc)

ols_step_best_subset(biochem_CHL_model_mc)

biochem_CHL_model2_mc <- lm(mag ~ L_CHL + E_CHL + A_CHL + 0, data = biochem_fit.coef2)
summary(biochem_CHL_model2_mc)

qqnorm(biochem_CHL_model2_mc$residuals)



#Adding (standard) Deviations into Quality Predictors model

Biochem_comp2_dev <- read_csv("Biochem_comp2_dev.csv")














################################# Morfi (Biomass / DL) Covariate ###########################


#Testing new Quality Predictor values - MorFi_DL and MorFi_Biomass

    ### MorFi (Biomass) correlation

morfi_biomass <- read_csv("morfi_biomass.csv")

morfi_biomass$biochem.ID <- paste(morfi_biomass$Variety, morfi_biomass$Section, sep = "-")

fits.coef$biochem.ID <- rownames(fits.coef)

fits.coef$mag <- sqrt((fits.coef$intercept^2)+(fits.coef$slope^2))

morfi_biomass2 <- merge(morfi_biomass,fits.coef,by="biochem.ID")

morfi_biomass2$mc <- morfi_biomass2$slope + morfi_biomass2$intercept



morfi_biomass_model <- lm(morfi_biomass2$mc ~ . + 0 , data = morfi_biomass2[,c(8:17,19:22,26:27)])

ols_step_all_possible(morfi_biomass_model)

ols_step_best_subset(morfi_biomass_model)

morfi_biomass_model2 <- lm(mc ~ fibre_L + fibre_broken + fine_cont.L.L + fine_A + 0, data = morfi_biomass2)
summary(morfi_biomass_model2)


    ### MorFi (DL) correlation

morfi_dl <- read_csv("morfi_dl.csv")

morfi_dl$biochem.ID <- paste(morfi_dl$Variety, morfi_dl$Section, sep = "-")

morfi_dl2 <- merge(morfi_dl,fits.coef,by="biochem.ID")

morfi_dl2$mc <- morfi_dl2$slope + morfi_dl2$intercept



morfi_dl_model <- lm(morfi_dl2$mc ~ . + 0 , data = morfi_dl2[,c(8:17,19:22,26:27)])

ols_step_all_possible(morfi_dl_model)

ols_step_best_subset(morfi_dl_model)

morfi_dl_model2 <- lm(mc ~ fibre_cont + fibre_L + fibre_L.L + fibre_L3.L + fibre_A.L + fibre_W + fibre_course + fibre_kink.num + fibre_kink.ang + MF.index + fibre_broken + fine_n + fine_A + 0, data = morfi_biomass2)
summary(morfi_dl_model2)



#Get average values from morfi_biomass

morfi_biomass_ave <- morfi_biomass %>% group_by(biochem.ID) %>% 
  summarize_at(.vars = names(.)[6:26], .funs = c(mean="mean"))

write.csv(morfi_biomass_ave,"C:\\Users\\s4294173\\Desktop\\Sorghum2020\\morfi_biomass_ave.csv", row.names = FALSE)











################################# Biochemical Comp Facets ###########################


biochem_comp_symbols <- read_csv("biochem comp symbols.csv")

biochem_comp_symbols$Biochemical = factor(biochem_comp_symbols$Biochemical, levels=c("Cellulose", "Hemicellulose","Lignin","Ash","Extractives"), labels=c("Cellulose", "Hemicellulose","Lignin","Ash","Extractives"))


p <- ggplot(biochem_comp_symbols) +
  aes(x=Variety, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(15,17,16,1,10)) +
  scale_color_manual(values=c("darkgreen", "orange2", "brown4","grey50","pink")) +
  geom_segment(aes(x=Variety, 
                   xend=Variety, 
                   y=min(Content), 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Biomass Biochemical Composition", 
       subtitle="Breakdown by sorghum variety and plant section") +  
  coord_flip()

p + facet_grid(rows = vars(Biochemical)) +
  theme_bw()












################################# TGA Covariate ###########################


    ### Correlation - TGA vs Nanopaper

TGA_sorg <- read_csv("TGA_sorghum.csv")

TGA_sorg$biochem.ID <- paste(TGA_sorg$Variety, TGA_sorg$Section, sep = "-")

TGA_sorg.coef <- merge(TGA_sorg,fits.coef,by="biochem.ID")

TGA_sorg.coef$mc <- TGA_sorg.coef$slope + TGA_sorg.coef$intercept


TGA_model_mc <- lm(mc ~ Step + Onset + Peak1 + Peak2 + Residue + 0, data = TGA_sorg.coef)

ols_step_all_possible(TGA_model_mc)

ols_step_best_subset(TGA_model_mc)

TGA_model2_mc <- lm(mc ~ Step + Onset + Peak1 + Residue + 0, data = TGA_sorg.coef)
summary(TGA_model2_mc)

qqnorm(TGA_model2_mc$residuals)













################################# XRD Covariate ###########################


### Correlation - XRD vs Nanopaper

XRD_sorg <- read_csv("XRD_sorghum.csv")

XRD_sorg$biochem.ID <- paste(XRD_sorg$Variety, XRD_sorg$Section, sep = "-")

XRD_sorg.coef <- merge(XRD_sorg,fits.coef,by="biochem.ID")

XRD_sorg.coef$mc <- XRD_sorg.coef$slope + XRD_sorg.coef$intercept


XRD_model_mc <- lm(mc ~ CrI + I_200 + I_amp + 0, data = XRD_sorg.coef)

ols_step_all_possible(XRD_model_mc)

ols_step_best_subset(XRD_model_mc)

XRD_model2_mc <- lm(mc ~ CrI + 0, data = XRD_sorg.coef)

summary(XRD_model2_mc)

qqnorm(XRD_model2_mc$residuals)











################################# Biomass Benchmarking BIOCHEM COMP ###########################

## Analysing covariate data - Biochemical composition

Biochem_BiomassBenchmarking2 <- read_csv("Biochem_BiomassBenchmarking2.csv")


## Generating new variable based on wt% basis and CHL basis

Biochem_BiomassBenchmarking2$C_wt <- Biochem_BiomassBenchmarking2$Glucan
Biochem_BiomassBenchmarking2$H_wt <- Biochem_BiomassBenchmarking2$Xylan + Biochem_BiomassBenchmarking2$Mannan + Biochem_BiomassBenchmarking2$Arabinan + Biochem_BiomassBenchmarking2$Galactan + Biochem_BiomassBenchmarking2$Rhamnan
Biochem_BiomassBenchmarking2$L_wt <- Biochem_BiomassBenchmarking2$Klason_Lignin + Biochem_BiomassBenchmarking2$AS_Lignin
Biochem_BiomassBenchmarking2$E_wt <- Biochem_BiomassBenchmarking2$Extractives
Biochem_BiomassBenchmarking2$A_wt <- Biochem_BiomassBenchmarking2$Ash + Biochem_BiomassBenchmarking2$Acid_Insoluble_Ash

Biochem_BiomassBenchmarking2$C_CHL <- Biochem_BiomassBenchmarking2$C_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt)
Biochem_BiomassBenchmarking2$H_CHL <- Biochem_BiomassBenchmarking2$H_wt/(Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt)
Biochem_BiomassBenchmarking2$L_CHL <- Biochem_BiomassBenchmarking2$L_wt/(Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt)

Biochem_BiomassBenchmarking2$C_CHLEA <- Biochem_BiomassBenchmarking2$C_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$E_wt + Biochem_BiomassBenchmarking2$A_wt)
Biochem_BiomassBenchmarking2$H_CHLEA <- Biochem_BiomassBenchmarking2$H_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$E_wt + Biochem_BiomassBenchmarking2$A_wt)
Biochem_BiomassBenchmarking2$L_CHLEA <- Biochem_BiomassBenchmarking2$L_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$E_wt + Biochem_BiomassBenchmarking2$A_wt)
Biochem_BiomassBenchmarking2$E_CHLEA <- Biochem_BiomassBenchmarking2$E_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$E_wt + Biochem_BiomassBenchmarking2$A_wt)
Biochem_BiomassBenchmarking2$A_CHLEA <- Biochem_BiomassBenchmarking2$A_wt/(Biochem_BiomassBenchmarking2$C_wt + Biochem_BiomassBenchmarking2$H_wt + Biochem_BiomassBenchmarking2$L_wt + Biochem_BiomassBenchmarking2$E_wt + Biochem_BiomassBenchmarking2$A_wt)

Biochem_BiomassBenchmarking2$C.H <- Biochem_BiomassBenchmarking2$C_CHLEA/Biochem_BiomassBenchmarking2$H_CHLEA
Biochem_BiomassBenchmarking2$C.L <- Biochem_BiomassBenchmarking2$C_CHLEA/Biochem_BiomassBenchmarking2$L_CHLEA
Biochem_BiomassBenchmarking2$C.E <- Biochem_BiomassBenchmarking2$C_CHLEA/Biochem_BiomassBenchmarking2$E_CHLEA
Biochem_BiomassBenchmarking2$C.A <- Biochem_BiomassBenchmarking2$C_CHLEA/Biochem_BiomassBenchmarking2$A_CHLEA
Biochem_BiomassBenchmarking2$H.L <- Biochem_BiomassBenchmarking2$H_CHLEA/Biochem_BiomassBenchmarking2$L_CHLEA
Biochem_BiomassBenchmarking2$H.E <- Biochem_BiomassBenchmarking2$H_CHLEA/Biochem_BiomassBenchmarking2$E_CHLEA
Biochem_BiomassBenchmarking2$H.A <- Biochem_BiomassBenchmarking2$H_CHLEA/Biochem_BiomassBenchmarking2$A_CHLEA
Biochem_BiomassBenchmarking2$L.E <- Biochem_BiomassBenchmarking2$L_CHLEA/Biochem_BiomassBenchmarking2$E_CHLEA
Biochem_BiomassBenchmarking2$L.A <- Biochem_BiomassBenchmarking2$L_CHLEA/Biochem_BiomassBenchmarking2$A_CHLEA
Biochem_BiomassBenchmarking2$E.A <- Biochem_BiomassBenchmarking2$E_CHLEA/Biochem_BiomassBenchmarking2$A_CHLEA





Biochem_BiomassBenchmarking2$biochem.ID <- paste(Biochem_BiomassBenchmarking2$Biomass, 
                                                 Biochem_BiomassBenchmarking2$Section,
                                                 Biochem_BiomassBenchmarking2$Treatment,
                                                 sep = "-")

Biochem_BiomassBenchmarking2$ID <- paste(Biochem_BiomassBenchmarking2$Biomass, 
                                                 Biochem_BiomassBenchmarking2$Section,
                                                 sep = "-")


# write.csv(Biochem_BiomassBenchmarking2,'biochem_biomass_symbols.csv') 
# 
# Biochem_Biomass_Symbols <- read_csv("Biochem_Biomass_Symbols.csv")


Biochem_Biomass_Symbols <- pivot_longer(Biochem_BiomassBenchmarking2,
                                        C_wt:E.A,
                                        names_to = "Biochemical",
                                        values_to = "Content")


Biochem_Biomass_Symbols$Biochemical = factor(Biochem_Biomass_Symbols$Biochemical,
                                             levels=c("C_wt", "H_wt","L_wt","E_wt","A_wt"),
                                             labels=c("Cellulose", "Hemicellulose","Lignin","Extractives","Ash"))


Biochem_Biomass_Symbols$Treatment = factor(Biochem_Biomass_Symbols$Treatment, 
                                             levels=c("Raw","DL"), 
                                             labels=c("Raw", "DL"))



p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Cellulose", "Hemicellulose","Lignin", "Ash"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("darkgreen", "orange2", "brown4","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Biochemical Composition", 
       subtitle="Biomass Benchmarking Study - Comparison across Biochemical Components, Biomass Types, Sections, and Treatment Types") +  
  coord_flip()

p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw()






p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Cellulose", "Hemicellulose","Lignin", "Ash"))) +
  aes(x=ID, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("darkgreen", "orange2", "brown4","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Biochemical Composition", 
       subtitle="Biomass Benchmarking Study - Comparison across Biochemical Components, Biomass Types, Sections, and Treatment Types") +  
  coord_flip()

p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw()




p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Extractives","Ash"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("bisque3","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Biochemical Composition", 
       subtitle="Biomass Benchmarking Study") +  
  coord_flip()

p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw()





# p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("C.H", "C.L","H.L"))) +
#   aes(x=Biomass, y=Content) + 
#   geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
#   scale_shape_manual(values=c(12,
#                               15,
#                               7,
#                               17,
#                               19,
#                               1,
#                               10)) +
#   scale_color_manual(values=c("darkgreen", "orange2", "brown4","bisque3","grey50")) +
#   geom_segment(aes(x=Biomass, 
#                    xend=Biomass, 
#                    y=0, 
#                    yend=max(Content)), 
#                linetype="dashed", 
#                size=0.1) +   # Draw dashed lines
#   labs(title="Biochemical Composition", 
#        subtitle="Biomass Benchmarking Study") +  
#   coord_flip()
# 
# p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
#   theme_bw()



########################################################################################

### Cellulose

p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Cellulose"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("darkgreen", "orange2", "brown4","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title="Biochemical Composition", 
       subtitle="Biomass Benchmarking Study") +  
  coord_flip()

C <- p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw() +
  theme(strip.text.y = element_blank()) + 
  guides(color = "none", shape = "none")


### Hemicellulose

p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Hemicellulose"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("orange2", "brown4","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title=" ", 
       subtitle=" ") +  
  coord_flip()

H <- p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw() + 
  theme(strip.text.y = element_blank(), 
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
  guides(color = FALSE, shape = FALSE)




### Lignin

p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Lignin"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("brown4","grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title=" ", 
       subtitle=" ") +  
  coord_flip()

L <- p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw() + 
  theme(strip.text.y = element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(color = FALSE, shape = FALSE)


### Ash

p <- ggplot(subset(Biochem_Biomass_Symbols, Biochemical %in% c("Ash"))) +
  aes(x=Biomass, y=Content) + 
  geom_point(size=5, aes(shape = Section, color = Biochemical)) +   # Draw points
  scale_shape_manual(values=c(12,
                              15,
                              7,
                              17,
                              19,
                              1,
                              10)) +
  scale_color_manual(values=c("grey50")) +
  geom_segment(aes(x=Biomass, 
                   xend=Biomass, 
                   y=0, 
                   yend=max(Content)), 
               linetype="dashed", 
               size=0.1) +   # Draw dashed lines
  labs(title=" ", 
       subtitle=" ") +  
  coord_flip()

A <- p + facet_grid(rows = vars(Treatment), cols = vars(Biochemical)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(color = FALSE)



grid.arrange(C, H, L, A, nrow=1) 







