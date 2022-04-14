

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(lme4)
library(emmeans)
library(effects)
library(rstatix)
library(lmerTest)
library(selection.index)
library(vegan)
library(ggpubr)
library(scales)
library(ggthemes)
library(gghighlight)



############ .csv Files #################


sedimentation <- read_csv("sedimentation.csv")

WRV <- read_csv("WRV.csv")

morfi_biomass <- read_csv("morfi_biomass.csv")

morfi_cnf2<- read_csv("morfi_cnf2.csv")


colnames(morfi_Q)[10] <- "fibre_coarse"



############ Regression Slope Matrix #################

morfi_Q_long <- gather(morfi_Q, parameter, value, fibre_n:fine_L, factor_key=TRUE)


morfi_quality_joined <- right_join(morfi_Q_long,qr_scaling2, by="ID")

fits.morfi_quality <- lmList(score.y ~ as.numeric(value) | with(morfi_quality_joined, interaction(quality_def,parameter)), data = morfi_quality_joined)

fits.coef.morfi_quality <- coef(fits.morfi_quality)
colnames(fits.coef.morfi_quality) <- c("intercept","slope")

fits.coef.morfi_quality$ID <- rownames(fits.coef.morfi_quality)

ID2 <- str_split(fits.coef.morfi_quality$ID,"\\.",simplify = TRUE, n=2)

fits.coef.morfi_quality <- cbind(ID2,fits.coef.morfi_quality)

colnames(fits.coef.morfi_quality) <- c("quality_def","parameter","intercept","slope","ID")



### Q1 ###

fits.coef.morfi_quality.Q1 <- subset(fits.coef.morfi_quality, quality_def == "Q1")
fits.coef.morfi_quality.Q1$slope_direction <- ifelse(fits.coef.morfi_quality.Q1$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q1 <- fits.coef.morfi_quality.Q1[order(fits.coef.morfi_quality.Q1$slope), ]
fits.coef.morfi_quality.Q1$parameter <- factor(fits.coef.morfi_quality.Q1$parameter, levels = fits.coef.morfi_quality.Q1$parameter)



beta_Q1 <- ggplot(fits.coef.morfi_quality.Q1, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q1) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()




### Q2 ###

fits.coef.morfi_quality.Q2 <- subset(fits.coef.morfi_quality, quality_def == "Q2")
fits.coef.morfi_quality.Q2$slope_direction <- ifelse(fits.coef.morfi_quality.Q2$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q2 <- fits.coef.morfi_quality.Q2[order(fits.coef.morfi_quality.Q2$slope), ]
fits.coef.morfi_quality.Q2$parameter <- factor(fits.coef.morfi_quality.Q2$parameter, levels = fits.coef.morfi_quality.Q2$parameter)



beta_Q2 <- ggplot(fits.coef.morfi_quality.Q2, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q2) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()




### Q3 ###

fits.coef.morfi_quality.Q3 <- subset(fits.coef.morfi_quality, quality_def == "Q3")
fits.coef.morfi_quality.Q3$slope_direction <- ifelse(fits.coef.morfi_quality.Q3$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q3 <- fits.coef.morfi_quality.Q3[order(fits.coef.morfi_quality.Q3$slope), ]
fits.coef.morfi_quality.Q3$parameter <- factor(fits.coef.morfi_quality.Q3$parameter, levels = fits.coef.morfi_quality.Q3$parameter)



beta_Q3 <- ggplot(fits.coef.morfi_quality.Q3, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q3) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()




### Q4 ###

fits.coef.morfi_quality.Q4 <- subset(fits.coef.morfi_quality, quality_def == "Q4")
fits.coef.morfi_quality.Q4$slope_direction <- ifelse(fits.coef.morfi_quality.Q4$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q4 <- fits.coef.morfi_quality.Q4[order(fits.coef.morfi_quality.Q4$slope), ]
fits.coef.morfi_quality.Q4$parameter <- factor(fits.coef.morfi_quality.Q4$parameter, levels = fits.coef.morfi_quality.Q4$parameter)



beta_Q4 <- ggplot(fits.coef.morfi_quality.Q4, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q4) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()




### Q5 ###

fits.coef.morfi_quality.Q5 <- subset(fits.coef.morfi_quality, quality_def == "Q5")
fits.coef.morfi_quality.Q5$slope_direction <- ifelse(fits.coef.morfi_quality.Q5$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q5 <- fits.coef.morfi_quality.Q5[order(fits.coef.morfi_quality.Q5$slope), ]
fits.coef.morfi_quality.Q5$parameter <- factor(fits.coef.morfi_quality.Q5$parameter, levels = fits.coef.morfi_quality.Q5$parameter)



beta_Q5 <- ggplot(fits.coef.morfi_quality.Q5, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q5) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()



### Q6 ###

fits.coef.morfi_quality.Q6 <- subset(fits.coef.morfi_quality, quality_def == "Q6")
fits.coef.morfi_quality.Q6$slope_direction <- ifelse(fits.coef.morfi_quality.Q6$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q6 <- fits.coef.morfi_quality.Q6[order(fits.coef.morfi_quality.Q6$slope), ]
fits.coef.morfi_quality.Q6$parameter <- factor(fits.coef.morfi_quality.Q6$parameter, levels = fits.coef.morfi_quality.Q6$parameter)



beta_Q6 <- ggplot(fits.coef.morfi_quality.Q6, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q6) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()



### Q7 ###

fits.coef.morfi_quality.Q7 <- subset(fits.coef.morfi_quality, quality_def == "Q7")
fits.coef.morfi_quality.Q7$slope_direction <- ifelse(fits.coef.morfi_quality.Q7$slope < 0, 
                                                     "negative", 
                                                     "positive")
fits.coef.morfi_quality.Q7 <- fits.coef.morfi_quality.Q7[order(fits.coef.morfi_quality.Q7$slope), ]
fits.coef.morfi_quality.Q7$parameter <- factor(fits.coef.morfi_quality.Q7$parameter, levels = fits.coef.morfi_quality.Q7$parameter)



beta_Q7 <- ggplot(fits.coef.morfi_quality.Q7, aes(x=parameter, y =slope, label=slope)) + 
  geom_bar(stat='identity', aes(fill=slope_direction), width=.5)  +
  xlab("MorFi Parameters") +
  ylab("Quality Gradient") +
  scale_fill_manual(name="Influence on Quality", 
                    labels = c("Positive", "Negative"), 
                    values = c("positive"="#00ba38", "negative"="#f8766d")) + 
  labs(subtitle="Quality Gradient (Q7) for All MorFi Parameters", 
       title= "Diverging Bars") + 
  facet_grid(~quality_def) +
  coord_flip() +
  ylim(-16,16) +
  theme_bw()


### Arrange all Quality Gradient Plots ###

ggarrange(beta_Q1, beta_Q2, beta_Q4, beta_Q5, beta_Q6, beta_Q7,
          # labels = c("A", "B", "C"),
          ncol = 3, nrow = 2,
          common.legend = TRUE)


test <- merge(x = fits.coef.morfi_quality.Q1, 
              y = fits.coef.morfi_quality.Q2, 
              by = "parameter", all.x = TRUE)

# test <- merge(x = test, 
#               y = fits.coef.morfi_quality.Q7, 
#               by = "parameter", all.x = TRUE)








# beta_Q7_lollypop <- ggplot(fits.coef.morfi_quality.Q7, aes(x=parameter, y =round(slope, digits = 0), label=slope)) +
#   geom_point(stat='identity', aes(color =slope_direction),  size=9)  +
#   # geom_segment(aes(y = 0,
#   #                  x = parameter,
#   #                  yend = slope,
#   #                  xend = parameter,
#   #                  color= slope_direction)) +
#   scale_color_manual(name="Influence on Quality",
#                      labels = c("Positive", "Negative"),
#                      values = c("positive"="#00ba38", "negative"="#f8766d")) +
#   geom_text(color="white", size=4) +
#   xlab("MorFi Parameter") +
#   ylab("Quality Gradient") +
#   labs(subtitle="Quality Gradient (Q7) for All MorFi Parameters",
#        title= "Diverging Bars") +
#   facet_grid(~quality_def) +
#   coord_flip() +
#   theme_bw()
# 
# beta_Q7_lollypop



fits.coef.morfi_quality1 <- fits.coef.morfi_quality %>% group_by(parameter) %>% dplyr::mutate(mean=mean(slope))


fits.coef.morfi_quality2 <- fits.coef.morfi_quality1 %>% dplyr::mutate(parameter_order = fct_reorder(parameter, mean,.desc=TRUE))


g <- ggplot(fits.coef.morfi_quality2, aes(x=fct_reorder(parameter_order,mean), y =slope, fill=parameter))
g + geom_violin(aes(fill = parameter), size = 1.25) + 
  labs(title="Violin plot - Quality", 
       subtitle="Fibre Morphology vs Quality Gradients (Q1-7)",
       x="MorFi Parameter",
       y="Quality Gradient") +
  coord_flip() + 
  # gghighlight(max(value) > 20)
  # facet_grid(~quality_def) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14))






ggplot(beta_Q1Q5_mean, aes(fct_reorder(parameter_order,mean), y =slope, fill=parameter)) +
  geom_boxplot() +
  geom_point(aes(shape = quality_def,colour=parameter), size = 3.5) +
  labs(title="Quality Gradient Boxplots", 
       subtitle="Comparison between Q1 vs Q5",
       x="MorFi Parameter",
       y="Quality Gradient") +
  coord_flip() + 
  # facet_grid(~quality_def) +
  theme_minimal()












beta_Q1Q5 <- dplyr::filter(fits.coef.morfi_quality,grepl("Q1|Q5",quality_def))

beta_Q1Q5_mean <- beta_Q1Q5 %>% dplyr::group_by(parameter) %>% dplyr::mutate(mean=mean(slope)) %>%
  dplyr::mutate(parameter_order = fct_reorder(parameter, mean,.fun='mean'))


ggplot(beta_Q1Q5_mean, aes(fct_reorder(parameter_order,mean), y =slope, fill=parameter)) +
  geom_boxplot() +
  geom_point(aes(shape = quality_def,colour=parameter), size = 3.5) +
  labs(title="Quality Gradient Boxplots", 
       subtitle="Comparison between Q1 vs Q5",
       x="MorFi Parameter",
       y="Quality Gradient") +
  coord_flip() + 
  # facet_grid(~quality_def) +
  theme_minimal()













