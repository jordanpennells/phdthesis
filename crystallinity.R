
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(olsrr)
library(scales)
library(csv)



xrd <- read_csv("crystallinity_data.csv")


# x <- cbind(xrd$angle,xrd$angle,xrd$angle,xrd$angle,xrd$angle)
# y <- cbind(xrd$wood,xrd$banana_leaf,xrd$banana_stem,xrd$sugarcane_bagasse,xrd$sugarcane_mulch)
# 
# matplot(x,y,type="p")


d <- melt(xrd, id.vars="angle")

ggplot(d, aes(angle,value, col=variable)) + 
  # geom_point() + 
  geom_line(size = 1)

x <- ggplot(d, aes(angle,value, col=variable)) + 
  geom_line(size = 1) + 
  facet_wrap(~variable) +
  labs(title="Full XRD Panel",subtitle="Biomass Benchmarking Study") +
  theme_bw()

x + theme(legend.position = "none")



################################################

blank <- grid.rect(gp=gpar(col="white"))
p_void <- ggplot() + theme_void()
p_void <- grid.arrange(p_void, nrow=1,ncol=1) 


p <- ggplot(d, aes(angle,value, col=variable)) + 
  geom_line(size = 1) + 
  facet_wrap(~variable) +
  theme_bw()


p1 <- ggplot(subset(d, variable %in% c("wood","banana_leaf","banana_stem","sugarcane_mulch","sugarcane_bagasse"))) + 
         aes(angle,value, col=variable) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, nrow=1) +
  xlab(" ") +
  theme_bw()


p2 <- ggplot(subset(d, variable %in% c("sugargraze_leaf",
                                       "sugargraze_sheath",
                                       "sugargraze_lowerstem",
                                       "sugargraze_upperstem"))) + 
  aes(angle,value, col=variable) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, nrow=1) +
  xlab(" ") +
  theme_bw()


# ps1 <- ggplot(subset(d, variable %in% c("sugargraze_leaf"))) + 
#   aes(angle,value, col=variable) + 
#   geom_line(size = 1) + 
#   facet_wrap(~variable, nrow=1) +
#   xlab(" ") +
#   theme_bw()
# 
# ps2 <- ggplot(subset(d, variable %in% c("sugargraze_sheath"))) + 
#   aes(angle,value, col=variable) + 
#   geom_line(size = 1) + 
#   facet_wrap(~variable, nrow=1) +
#   xlab(" ") +
#   theme_bw()
# 
# 
# p2 <- grid.arrange(ps1,ps2,p_void, nrow=1) 



p3 <- ggplot(subset(d, variable %in% c("yemen_leaf",
                                       "yemen_sheath",
                                       "yemen_lowerstem",
                                       "yemen_upperstem"))) + 
  aes(angle,value, col=variable) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, nrow=1) +
  xlab(" ") +
  theme_bw()


p4 <- ggplot(subset(d, variable %in% c("greenleafBMR_leaf",
                                       "greenleafBMR_sheath",
                                       "greenleafBMR_stem"))) + 
  aes(angle,value, col=variable) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, nrow=1, ncol=3) +
  xlab(" ") +
  theme_bw()


p5 <- ggplot(subset(d, variable %in% c("graingrass_leaf",
                                       "graingrass_sheath",
                                       "graingrass_stem"))) + 
  aes(angle,value, col=variable) + 
  geom_line(size = 1) + 
  facet_wrap(~variable, nrow=1) +
  xlab(" ") +
  theme_bw()


# d$variable2 <- factor(d$variable,levels=sort(c("wood","banana_leaf","banana_stem","sugarcane_mulch","sugarcane_bagasse",
#                                           "sugargraze_leaf","sugargraze_sheath","sugargraze_lowerstem","sugargraze_upperstem"," ",
#                                           "yemen_leaf",
#                                           "yemen_sheath",
#                                           "yemen_lowerstem",
#                                           "yemen_upperstem", " ",
#                                           "greenleafBMR_leaf",
#                                           "greenleafBMR_sheath",
#                                           "greenleafBMR_stem", " ", " ",
#                                           "graingrass_leaf",
#                                           "graingrass_sheath",
#                                           "graingrass_stem"," "," ")))
# 
# d$variable2 <- factor(d$variable,levels=c("","",unique(d$variable)))


ggplot(d, aes(angle,value, col=variable)) + 
  geom_line(size = 1) + 
  facet_wrap(~variable) +
  theme_bw()






grid.arrange(p1,p2,p_void,p3,p4,p5, nrow=5) 



# g = ggplotGrob(p)
# 
# g$grobs[names(g$grobs) %in% c("panel1", "panel2", "strip_t.1", "strip_t.2")] = NULL
# g$layout = g$layout[!(g$layout$name %in% c("panel-1", "panel-2", 
#                                            "strip_t-1", "strip_t-2")),]
# g$layout[g$layout$name == "axis_l-1", c("l", "r")] = c(9,9)
# grid.newpage()
# grid.draw(g)















