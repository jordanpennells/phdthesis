

##Categorical scatterplot for Quality Ranking data

library(readr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggfittext)
library(ggrepel)
theme_set(theme_pubr())



qr <- read_csv("quality_ranking.csv")

qr$HPH.fac = factor(qr$HPH, levels=c("L", "M", "H"))

qr$ID = str_c(qr$variety,'-',qr$section,'-',qr$HPH)

#ggplot(qr, aes(x = quality_def, y = score)) + geom_bar(aes(fill = variety), stat = "identity", color = "white", position = position_dodge(0.9)) + facet_wrap(~HPH.fac) + fill_palette("jco")

#qr$HPH = as.numeric(levels(qr$HPH))[qr$HPH]
#qr$variety=as.numeric(levels(qr$variety))[qr$variety]
#ggballoonplot(qr, x = "quality_def", y = "score", size = "HPH", fill = "variety", ggtheme = theme_bw()) + scale_fill_viridis_c(option = "C") + scale_size_continuous(limits=c(400, 1100), breaks=c(400, 700, 1100)) + scale_fill_continuous(limits=c(1, 4), breaks=c(1, 2, 3, 4))


#ggballoonplot(qr, x = "quality_def", y = "score", shape = "section", size = "HPH.fac", fill = "variety", ggtheme = theme_bw()) + scale_size_discrete() + scale_fill_discrete() 

#vars <- c('Q1', 'Q2', 'Q3')

#df <- qr[!is.na(qr$variety),c('section','variety',vars)] %>%
  #group_by(variety) %>% sample_n(5) %>% ungroup() %>%
  #gather(key=quality_def,value=score,-section,-variety) %>%
  #filter(!is.na(score)) %>% 
  #arrange (variety, section)

#vars_x_axis <- c(qr %>% arrange(col) %>% select(quality_def) %>% distinct())$quality_def


  #Quality Ranking graph
ggplot(qr) +
  aes(x = quality_def, y = score, color = variety, shape = section, size = HPH.fac) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point() +
  geom_text(aes(label=ID),hjust=1.2,vjust=0,size=3) +
  theme(legend.position = "right", 
        panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Quality Definition") +
  ylab("Quality Score")



  #Quality Ranking graph (top/bottom points labelled)


S <- as.data.frame(aggregate(score~quality_def,data = qr,FUN = quantile,0.9))

colnames(S) <- c("quality_def","ID3")

qr2 <- qr %>% 
  full_join(S, by = "quality_def") %>%
  group_by(quality_def) %>%
  mutate(ID2 = ifelse(score>=ID3,ID,"")) %>%
  ungroup()


ggplot(qr2) +
  aes(x = quality_def, y = score, color = variety, shape = section, size = HPH.fac) +
  labs(shape = "Plant Section\n", color = "Sorghum Variety\n", size = "HPH Energy Level\n") +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_label(aes(label=ID2),hjust=-0.05,vjust=0,size=3.8, fontface = "bold", position=position_jitterdodge(jitter.width=0,jitter.height=1, dodge.width = 0, seed = 1), alpha = 0,label.size = 0) +
  theme(axis.text.y = element_text(face="bold"),axis.text.x = element_text(face="bold"),legend.position = "right", legend.text=element_text(size=12),legend.title=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Quality Definition") +
  ylab("Quality Score") +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  expand_limits(x=c(0.7,6), y=c(20, 100))




Q1 <- qr %>%  filter(str_detect(quality_def,"Q1"))
Q2 <- qr %>%  filter(str_detect(quality_def,"Q2"))
Q3 <- qr %>%  filter(str_detect(quality_def,"Q3"))
Q4 <- qr %>%  filter(str_detect(quality_def,"Q4"))
Q5 <- qr %>%  filter(str_detect(quality_def,"Q5"))


ggplot(qr2) +
  aes(x = 1, y = score, color = variety, shape = section, size = HPH.fac) +
  labs(shape = "Plant Section\n", color = "Sorghum Variety\n", size = "HPH Energy Level\n") +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point() +
  facet_grid(~quality_def) +
  theme(strip.text.x = element_text(size = 14, face = "bold"),axis.text.y = element_text(size = 14, face="bold"),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "right", legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Quality Definition") +
  ylab("Quality Score") +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  expand_limits(x=c(0.7,6), y=c(20, 100)) +
  geom_text_repel(aes(label=ID2), size = 5.5, fontface = "bold", hjust = -0.5, direction = "y", min.segment.length = 0)








##################### Scaling Factor #####################


qr_scaling <- read_csv("quality_scaling.csv")

qr_scaling$HPH.fac = factor(qr_scaling$HPH, levels=c("L", "M", "H"))

qr_scaling$ID = str_c(qr_scaling$variety,'-',qr_scaling$section,'-',qr_scaling$HPH)


S_scaling <- as.data.frame(aggregate(score~quality_def,data = qr_scaling,FUN = quantile,0.9))

colnames(S_scaling) <- c("quality_def","ID3")

qr_scaling2 <- qr_scaling %>% 
  full_join(S_scaling, by = "quality_def") %>%
  group_by(quality_def) %>%
  mutate(ID2 = ifelse(score>=ID3,ID,"")) %>%
  ungroup()

Q1 <- qr_scaling %>%  filter(str_detect(quality_def,"Q1"))
Q2 <- qr_scaling %>%  filter(str_detect(quality_def,"Q2"))
Q3 <- qr_scaling %>%  filter(str_detect(quality_def,"Q3"))
Q4 <- qr_scaling %>%  filter(str_detect(quality_def,"Q4"))
Q5 <- qr_scaling %>%  filter(str_detect(quality_def,"Q5"))
Q6 <- qr_scaling %>%  filter(str_detect(quality_def,"Q6"))
Q7 <- qr_scaling %>%  filter(str_detect(quality_def,"Q7"))


ggplot(subset(qr_scaling2, quality_def %in% c("Q1","Q5")), aes(x = 1, y = score, color = variety, shape = section, size = HPH.fac)) +
  labs(shape = "Plant Section\n", color = "Sorghum Variety\n", size = "HPH Energy Level\n") +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) +
  geom_point() +
  facet_grid(.~quality_def) +
  theme(strip.text.x = element_text(size = 14, face = "bold"),axis.text.y = element_text(size = 14, face="bold"),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "right", legend.text=element_text(size=16),legend.title=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"), panel.background = element_rect(fill = "white",color="black")) + 
  xlab("Quality Definition") +
  ylab("Quality Score") +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  expand_limits(x=c(0.7,6), y=c(20, 100)) +
  geom_text_repel(aes(label=ID2), size = 5.5, fontface = "bold", hjust = -0.5, direction = "y", min.segment.length = 0)






##################### Energy vs Quality #####################


qr_scaling2$HPH.number <- gsub("L",0,gsub("M",0.5,gsub("H",1,qr_scaling2$HPH)))


ggplot(qr_scaling2) +
  aes(x = HPH.number, y = score, color = variety) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange","darkgrey","red","yellow")) +
  # scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point(size=3.5, colour="black") + 
  stat_summary(aes(group=quality_def, colour=quality_def), fun=mean, geom="smooth") +
  ylab("Quality score")





ggplot(qr_scaling2) +
  aes(x = HPH.number, 
      y = score,  
      color = variety) +
  scale_color_manual(values=c("purple", "darkgreen",rainbow(7), "blue", "darkorange")) +
  # scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point(size=3.5) +
  stat_summary(aes(group=quality_def,
                   colour=quality_def), 
               fun=mean, 
               geom="smooth") +
  ylab("Quality score")



ggplot(qr_scaling2) +
  aes(x = HPH.number, 
      y = score,  
      color = variety) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  # scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point(size=3.5) +
  stat_summary(aes(group=quality_def,linetype = quality_def), 
               fun=mean,
               colour = "black",
               geom="smooth") +
  ylab("Quality score")



ggplot(qr_scaling2) +
  aes(x = HPH.number, 
      y = score,  
      color = variety, shape = section) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point(size=3.5) +
  stat_summary(aes(group=quality_def,linetype = quality_def), 
               fun=mean,
               colour = "black",
               geom="smooth") +
  ylab("Quality score")
















ggplot(qr_scaling2) +
  aes(x = HPH.number, 
      y = score,  
      color = variety,
      group=quality_def) +
  scale_color_manual(values=c("purple", "darkgreen", "blue", "darkorange")) +
  # scale_shape_manual(values=c(15,17,16,1,10)) + 
  theme_bw() + 
  geom_point(size=3.5,alpha = 0.2) +
  stat_summary(aes(linetype = quality_def), 
               fun=mean,
               colour = "black",
               geom="smooth",
               size=1) +
  #scale_linetype_identity() +
  scale_linetype(name = "Quality Definition")+
  ylab("Quality score") +
  theme(legend.text = element_text(colour="black", 
                                   size = 10, 
                                   face = "bold"), 
        legend.key.width = unit(1,"cm"))










fits.quality <- lmList(score ~ as.numeric(HPH.number) | quality_def, data = qr_scaling2)

fits.quality.coef <- coef(fits.quality)

colnames(fits.quality.coef) <- c("intercept","slope")

fits.quality.coef$Q <- rownames(fits.quality.coef)




































