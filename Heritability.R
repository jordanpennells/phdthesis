
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


################### Read CSV files ################### 

sedimentation <- read_csv("sedimentation.csv")

WRV <- read_csv("WRV.csv")

morfi_biomass <- read_csv("morfi_biomass.csv")

morfi_cnf2<- read_csv("morfi_cnf2.csv")


################### MorFi Parameter Heritability ###################


morfi_cnf2$ID = str_c(morfi_cnf2$Variety,'-',morfi_cnf2$Section)


# morfi_cnf_scaled <- as_tibble(scale(morfi_cnf2[,7:27], scale=TRUE,center = TRUE))


morfi_names <- colnames(morfi_cnf2[7:27])

heritability <- data.frame(name=character(),h2sv=numeric(),h2v=numeric(),h2hph=numeric(),h2=numeric())

 for (i in morfi_names){
  h.model <-  lmer( as.formula(paste(i,"~ (1|HPH) + (1|Variety/Section)")) , data=morfi_cnf2) 
  h.vc <- VarCorr(h.model)
  h.vc <- as.data.frame(h.vc)
  heritability <- heritability %>% add_row(name=i,
                           h2sv = h.vc$vcov[1]/sum(h.vc$vcov),
                           h2v = h.vc$vcov[2]/sum(h.vc$vcov),
                           h2hph = h.vc$vcov[3]/sum(h.vc$vcov),
                           h2 = (h.vc$vcov[1]+h.vc$vcov[2])/sum(h.vc$vcov))
  # h <- h + 1
 }





h2.round <- heritability %>% mutate(across(.cols = 2:5,round,digits=2))



h2.round.diff <- h2.round %>% mutate(diff_svv = h2sv - h2v, 
                                     diff_svhph = h2sv - h2hph, 
                                     diff_vhph = h2v - h2hph)





############ Heritability graphs #################

h2_long <- gather(h2.round, factor, heritability, h2sv:h2hph, factor_key=TRUE)

h2_long$factor.fac <- factor(h2_long$factor, levels = c("h2v", "h2hph","h2sv"))


ggdensity(h2_long, x = "heritability",
          add = "mean", rug = TRUE,
          color = "factor.fac", fill = "factor.fac",
          palette = c("#00AFBB", "#E7B800", "#90A2DA"))


gghistogram(h2_long, x = "heritability",
            add = "mean", rug = TRUE,
            color = "factor.fac", fill = "factor.fac",
            palette = c("#00AFBB", "#E7B800","#90A2DA"))





p <- ggboxplot(h2_long, x = "factor.fac", y = "heritability",
          color = "factor.fac", palette =c("#00AFBB", "#E7B800", "#90A2DA"),
          add = "jitter", shape = "factor.fac")


my_comparisons <- list( c("h2v", "h2hph"), c("h2hph", "h2sv"), c("h2v", "h2sv") )


p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1)  


ggviolin(h2_long, x = "factor.fac", y = "heritability", fill = "factor.fac",
         palette = c("#00AFBB", "#E7B800", "#90A2DA"), add = "mean", 
         add.params = list(size = 1.5))+
  # scale_y_continuous(limits = c(-0.5, 1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1)




ggbarplot(h2.round, x = "name", y = "h2sv",
          fill = "name",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # palette = c("#00AFBB", "#E7B800", "#90A2DA"),
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          title = "Heritability - section:variety",
          xlab = "MorFi parameter",
          ylab = "Heritability (h2)"
          )


ggbarplot(h2.round, x = "name", y = "h2v",
          fill = "name",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # palette = c("#00AFBB", "#E7B800", "#90A2DA"),
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          title = "Heritability - variety",
          xlab = "MorFi parameter",
          ylab = "Heritability (h2)")


ggbarplot(h2.round, x = "name", y = "h2hph",
          fill = "name",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # palette = c("#00AFBB", "#E7B800", "#90A2DA"),
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          title = "Heritability - HPH energy level",
          xlab = "MorFi parameter",
          ylab = "Heritability (h2)")


ggbarplot(h2.round, x = "name", y = "h2",
          fill = "name",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # palette = c("#00AFBB", "#E7B800", "#90A2DA"),
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          title = "Fibre Morphology Heritability",
          xlab = "MorFi Parameter",
          ylab = "Heritability score (h2)",
          legend = "NULL")





h2.round$h2sv_z <- (h2.round$h2sv -mean(h2.round$h2sv))/sd(h2.round$h2sv)
h2.round$h2v_z <- (h2.round$h2v -mean(h2.round$h2v))/sd(h2.round$h2v)
h2.round$h2hph_z <- (h2.round$h2hph -mean(h2.round$h2hph))/sd(h2.round$h2hph)

h2.round$h2sv_grp <- factor(ifelse(h2.round$h2sv_z < 0, "low", "high"), 
                      levels = c("low", "high"))
h2.round$h2v_grp <- factor(ifelse(h2.round$h2v_z < 0, "low", "high"), 
                            levels = c("low", "high"))
h2.round$h2hph_grp <- factor(ifelse(h2.round$h2hph_z < 0, "low", "high"), 
                            levels = c("low", "high"))

ggbarplot(h2.round, x = "name", y = "h2sv_z",
          fill = "h2sv_grp",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 0,          # Rotate vertically x axis texts
          ylab = "z-score - section:variety heritability ",
          xlab = "MorFi parameter",
          legend.title = "section:variety Group",
          title = "Deviation Graph - section:variety ",
          rotate = TRUE
)


ggbarplot(h2.round, x = "name", y = "h2v_z",
          fill = "h2v_grp",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 0,          # Rotate vertically x axis texts
          ylab = "z-score - variety heritability",
          xlab = "MorFi parameter",
          legend.title = "variety Group",
          title = "Deviation Graph - variety",
          rotate = TRUE
)


ggbarplot(h2.round, x = "name", y = "h2hph_z",
          fill = "h2hph_grp",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 0,          # Rotate vertically x axis texts
          ylab = "z-score - HPH heritability",
          xlab = "MorFi parameter",
          legend.title = "HPH Group",
          title = "Deviation Graph - HPH energy level",
          rotate = TRUE
)

# ggbarplot(h2_long, x = "factor.fac", y = "heritability",
#           color = "factor.fac", fill = "factor.fac",
#           palette = c("#00AFBB", "#E7B800","#90A2DA"))












############ Split by variety #################

heritability <- data.frame(name=character(),h2sv=numeric(),h2v=numeric(),h2hph=numeric())

for (i in morfi_names){
  h.model <-  lmer( as.formula(paste(i,"~ (1|HPH) + (1|Variety/Section)")) , data=morfi_cnf2) 
  h.vc <- VarCorr(h.model)
  h.vc <- as.data.frame(h.vc)
  # heritability$name[h] <- i
  # heritability$h2sv[h] <- h.vc$vcov[1]/sum(h.vc$vcov)
  # heritability$h2v[h] <- h.vc$vcov[2]/sum(h.vc$vcov)
  # heritability$h2hph[h] <- h.vc$vcov[3]/sum(h.vc$vcov)
  heritability <- heritability %>% add_row(name=i,
                                           h2sv = h.vc$vcov[1]/sum(h.vc$vcov),
                                           h2v = h.vc$vcov[2]/sum(h.vc$vcov),
                                           h2hph = h.vc$vcov[3]/sum(h.vc$vcov))
  # h <- h + 1
}







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








beta_Q7_lollypop <- ggplot(fits.coef.morfi_quality.Q7, aes(x=parameter, y =round(slope, digits = 0), label=slope)) +
  geom_point(stat='identity', aes(color =slope_direction),  size=9)  +
  # geom_segment(aes(y = 0,
  #                  x = parameter,
  #                  yend = slope,
  #                  xend = parameter,
  #                  color= slope_direction)) +
  scale_color_manual(name="Influence on Quality",
                    labels = c("Positive", "Negative"),
                    values = c("positive"="#00ba38", "negative"="#f8766d")) +
  geom_text(color="white", size=4) +
  xlab("MorFi Parameter") +
  ylab("Quality Gradient") +
  labs(subtitle="Quality Gradient (Q7) for All MorFi Parameters",
       title= "Diverging Bars") +
  facet_grid(~quality_def) +
  coord_flip() +
  theme_bw()

beta_Q7_lollypop


g <- ggplot(fits.coef.morfi_quality, aes(x=parameter, y =slope, fill=parameter))
g + geom_boxplot() + 
  labs(title="Violin plot", 
       subtitle="MorFi Parameter vs Quality Gradient",
       x="MorFi Parameter",
       y="Quality Gradient") +
  coord_flip() + 
  # gghighlight(max(value) > 20)
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




# grid.arrange(length.length_B,
#              length.length_L,
#              length.length_M,
#              length.length_H, 
#              ncol=1) 





# beta_Q1Q5 <- merge(fits.coef.morfi_quality.Q1,fits.coef.morfi_quality.Q5,by="parameter")
# 
# 
# brks <- seq(-15000000, 15000000, 5000000)
# lbls = paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")
# 
# 
# ggplot(beta_Q1Q5, aes(x=parameter, y =slope, fill = quality_def.x)) +   # Fill column
#   geom_bar(stat = "identity", width = .6) +   # draw the bars
#   scale_y_continuous(breaks = brks,   # Breaks
#                      labels = lbls) + # Labels
#   coord_flip() +  # Flip axes
#   labs(title="Email Campaign Funnel") +
#   theme_tufte() +  # Tufte theme from ggfortify
#   theme(plot.title = element_text(hjust = .5), 
#         axis.ticks = element_blank()) +   # Centre plot title
#   scale_fill_brewer(palette = "Dark2")







# ################### Read CSV files ################### 
# 
# genMat<- gen.varcov(data = morfi_cnf2[,7:27], genotypes = seldata[,2],
#                     replication = seldata[,6])







