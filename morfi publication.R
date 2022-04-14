

### Script for Morfi Publication

library(corrplot)
library(RColorBrewer)
library(kableExtra)
library(plyr)
library(dplyr)
library(seriation)
library(gplots)
library(tidyr)
library(tidyverse)
library(dendextend)
library(ggcorrplot)
library(gplots)
library(ggrepel)
library(qgraph)
library(diffcor)
library(psych)


# cor(morfi_cnf$fibre_n,morfi_cnf$fibre_L,use="complete.obs")
# plot(morfi_cnf$fibre_n,morfi_cnf$fibre_L)

pairs(morfi_cnf[,7:12])
pairs(morfi_cnf[,13:20])
pairs(morfi_cnf[,21:27])

pairs(morfi_cnf[,7:27])

pairs(nanopaper[,8:22])


colnames(morfi_cnf_PCA)[8] <- "fibre_coarse"

M = cor(morfi_cnf_PCA,use="complete.obs")

#colnames(M) <- c("FC",":K^o",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".")


# corrplot(M, order = 'FPC', type = 'lower')
# 
# png("corrplot.jpg", width = 680,height = 980)
# set.seed(123)


###################### Basic Corrplot ##########################


corrplot(M, order = 'FPC', type = 'lower',tl.cex = 1.2, diag = TRUE, col=brewer.pal(8,"RdYlGn"), tl.col="black", mar=c(1,1,1,1))

# dev.off()



###################### Qgrpah  ##########################

qgraph(M, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=5)

qgraph(M, shape="circle", posCol="darkgreen", negCol="darkred", layout="spring", vsize=5)



###################### Correlation value  ##########################

# 
# corCi(M, keys = NULL, n.iter = 100,  p = 0.05,overlap = FALSE, 
#       poly = FALSE, method = "pearson", plot=TRUE,minlength=5,n=NULL)
# 
# 
# 
# FisherZ <- function( r ) {
#   # PRECONDITIONS
#   if(missing(r)) return(NA)
#   if(any(r > 1) || any(r < -1)) 
#     stop("Fisher Z transformation is only defined for correlations (i.e. -1 <= r <= 1)")
#   
#   # Identify +/-1 because transformation cannot handle that
#   rp1 <- which(r == 1)
#   rm1 <- which(r == -1)
#   # Transform
#   z <- log( (1+r) / (1-r) ) / 2
#   # Set r= 1/-1 to Inf/-Inf
#   z[rp1] <- Inf
#   z[rm1] <- -Inf
#   #
#   return( z )
# }
# 
# FisherZInv <- function( z ) {
#   # PRECONDITIONS
#   if(missing(z)) return(NA)
#   if(any(is.nan(z))) z[which(is.nan(z))] <- 1 # 'FisherZ()' not defined for 1
#   
#   zp1 <- which(z == Inf) # z == plus 1
#   zm1 <- which(z ==-Inf) # z == minus 1
#   # Tranform inversely
#   r <- (exp(2*z)-1) / (exp(2*z)+1)
#   # Set r = 1/-1 to Inf/-Inf
#   r[zp1] <- 1
#   r[zm1] <- -1
#   return(r)
# }
# 
# 
# ZM <- FisherZ(M)
# MZ <- FisherZInv(ZM)
# round(ZM, 2)
# 
# ci <- corCi(Thurstone[1:6,1:6],n=213)
# 
# corCi(M, keys = NULL, n.iter = 100,  p = 0.05,overlap = FALSE,
#      poly = FALSE, method = "pearson", plot=TRUE,minlength=5,n=NULL)
# 
# corCi(M, n, conf.level = 0.95, alternative = c("two.sided", "less", "greater"))
# 
# n <- 30
# r <- seq(0, .9, .1)
# rc <- t(sapply(r, corCi, n=n))
# t <- r * sqrt(n-2) / sqrt(1-r^2)
# p <- (1 - pt(t, n-2)) / 2
# r.rc <- data.frame(r=r, z=FisherZ(r), lower=rc[,2], upper=rc[,3], t=t, p=p)
# 
# round(r.rc,2)






###################### Clustered Corrplot ##########################

corrplot(M, method = 'square', diag = FALSE, order = 'hclust', col=brewer.pal(8,"RdYlGn"), addrect = 3, tl.col="black", rect.col = 'black', rect.lwd = 3, tl.pos = 'd')




corrplot(M, 
         method = 'color', 
         addCoef.col = 'white', 
         diag = TRUE, 
         order = 'hclust', 
         col=colorRampPalette(c("#a50026","#f46d43", "white","white","white","white","#66bd63","#006837"))(8), 
         addrect = 4, 
         tl.col="black",
         tl.cex = 1.2,
         rect.col = 'black', 
         rect.lwd = 5, 
         tl.pos = 'bt',
         number.cex=0.7,
         cl.pos = "b",
         cl.cex = 1.2,
         cl.ratio = 0.1,
         mar=c(0,0,1,0),
         insig='blank',
         addgrid.col = NA)








# corrplot(M, method = 'color',p.mat = cor1[[1]], insig = "blank", addCoef.col="grey", 
# order = "AOE", tl.cex = 1/par("cex"),
# cl.cex = 1/par("cex"), addCoefasPercent = TRUE)
# par(cex = cex.before)




###################### Unique Corrplots ##########################

dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriate(d_corr, method = method, ...)
  i = get_order(s)
  return(i)
}

# Fast Optimal Leaf Ordering for Hierarchical Clustering
i = dist2order(M, 'OLO')
corrplot(M[i, i], cl.pos = 'n', col=brewer.pal(8,"RdYlGn"), tl.col="black")

# Quadratic Assignment Problem
i = dist2order(M, 'QAP_2SUM')
corrplot(M[i, i], cl.pos = 'n', col=brewer.pal(8,"RdYlGn"), tl.col="black")

### Multidimensional Scaling
i = dist2order(M, 'MDS_nonmetric')
corrplot(M[i, i], cl.pos = 'n',col=brewer.pal(8,"RdYlGn"), tl.col="black")

############# Simulated annealing
i = dist2order(M, 'ARSA')
corrplot(M[i, i], cl.pos = 'n', col=brewer.pal(8,"RdYlGn"), tl.col="black")

# # TSP solver
# i = dist2order(M, 'TSP')
# corrplot(M[i, i], cl.pos = 'n',col=brewer.pal(8,"RdYlGn"), tl.col="black")

# Spectral seriation
i = dist2order(M, 'Spectral')
corrplot(M[i, i], cl.pos = 'n',col=brewer.pal(8,"RdYlGn"), tl.col="black")


###################### Aesthetic Corrplot ##########################

testRes = cor.mtest(morfi_cnf_PCA, conf.level = 0.95)


# corrplot(M,
#          p.mat = testRes$p,
#          method = "color",
#          type = 'lower',
#          insig='blank',
#          addCoef.col="black",
#          order = "hclust",
#          hclust.method = "mcquitty",
#          number.cex=0.75,
#          tl.cex = 1.1,
#          col=brewer.pal(8,"RdYlGn"),
#          tl.col="black")

# corrplot(M,
#          p.mat = testRes$p,
#          method = "color",
#          insig='blank',
#          addCoef.col="black", 
#          order = "FPC", 
#          number.cex=0.75,
#          tl.cex = 1.1,
#          col=brewer.pal(4,"RdYlGn"), 
#          tl.col="black")


corr_alpha <- c("black", "darkgrey", "grey", "white", "white", "grey", "darkgrey", "black")

corrplot(M,
         p.mat = testRes$p,
         method = "color",
         type = 'lower',
         diag = TRUE,
         insig='blank',
         addCoef.col = "black", 
         order = "FPC", 
         number.cex=0.8,
         tl.cex = 1.1,
         col=brewer.pal(11,"RdYlGn"), 
         tl.col="black",
         cl.pos = 'b',
         mar=c(1,1,1,1))


corrplot(M,
         p.mat = testRes$p,
         method = "color",
         type = 'lower',
         diag = TRUE,
         insig='blank',
         addCoef.col = "white", 
         order = "FPC", 
         number.cex=0.8,
         tl.cex = 1.1,
         col=brewer.pal(8,"RdYlGn"), 
         tl.col="black",
         cl.pos = 'b',
         mar=c(1,1,1,1))

par(xpd=TRUE)

corrplot(M,
         p.mat = testRes$p,
         method = "color",
         type = 'lower',
         diag = TRUE,
         insig='blank',
         addCoef.col = "white",
         order = "FPC",
         number.cex=0.75,
         tl.cex = 1.1,
         col=colorRampPalette(c("#a50026","#f46d43", "#fee08b","white","white","#a6d96a","#66bd63","#006837"))(8),
         tl.col="black",
         cl.pos = 'b',
         mar=c(0,0,5,0))



corrplot(M,
         p.mat = testRes$p,
         method = "color",
         type = 'lower',
         diag = TRUE,
         insig='blank',
         addCoef.col = "white",
         order = "FPC",
         number.cex=0.8,
         tl.cex = 1.1,
         col=colorRampPalette(c("#a50026","#f46d43", "white","white","white","white","#66bd63","#006837"))(8),
         tl.col="black",
         cl.pos = 'b',
         mar=c(1,1,1,1))




###################### ggCorrPlot ##########################

ggcorrplot(M,
           p.mat = testRes$p,
           type = "upper",
           show.diag = TRUE,
           lab = TRUE,
           lab_alpha = TRUE,
           colors = c("red2", "white", "forestgreen"),
           outline.color = "black",
           insig = c("blank"),
           sig.level = 0.05)






###################### HEATMAP ##########################

morfi_cnf_heatmap <- read.csv("morfi_cnf_heatmap.csv")

morfi_cnf_heatmap <- morfi_cnf_heatmap %>% dplyr::group_by(Variety,HPH) %>% dplyr::summarise(across(.cols=everything(),mean,na.rm=TRUE))

morfi_cnf_heatmap2 <- morfi_cnf_heatmap[,3:23]

morfi_cnf_heatmap2 <- as_tibble(scale(morfi_cnf_heatmap2, scale=TRUE,center = TRUE))

morfi_cnf_heatmap <- bind_cols(morfi_cnf_heatmap[,1:2],morfi_cnf_heatmap2)

morfi_cnf_heatmap_num <- morfi_cnf_heatmap[,3:23]


heatmap.cor <- cor(morfi_cnf_heatmap_num, use="pairwise.complete.obs")

# typeof(heatmap.cor)



heatmap.dist <- as.dist(1 - heatmap.cor)
heatmap.tree <- hclust(heatmap.dist, method="complete")
plot(heatmap.tree)

heatmap.dend <- as.dendrogram(heatmap.tree)

clusters <- cutree(heatmap.dend, k=4, order_clusters_as_data = FALSE)
table(clusters)

clusters.df <- data.frame(ID = names(clusters), cluster = clusters)

cluster3.ID <- filter(clusters.df, cluster == 3)$ID

cat(as.character(cluster3.ID[1:5]), quote=FALSE,sep="\n")




heatmap.long <- pivot_longer(morfi_cnf_heatmap,cols=3:23)

heatmap.long <- heatmap.long %>% dplyr::group_by(Variety,HPH,name) %>% dplyr::summarise(value=mean(value,na.rm=TRUE))



color.scheme <- rev(brewer.pal(4,"RdYlGn"))

heatmap.long %>%
  ggplot(aes(x = HPH, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors=color.scheme, limits = c(-2.3,2.3)) + 
  xlab("Processing Energy Level ") +
  ylab("MorFi Parameter") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 18, vjust = 2),
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_text(size = 14, face="bold")) +
  guides(fill=guide_legend(title="  SD\nValue"))


heatmap.long %>%
  ggplot(aes(x = HPH, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours=colorRampPalette(c("#a50026","#f46d43", "#fee08b","white","white","#a6d96a","#66bd63","#006837"))(3), limits = c(-2.5,2.5)) + 
  xlab("Processing Energy Level ") +
  ylab("MorFi Parameter") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 18, vjust = 2),
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_text(size = 14, face="bold")) +
  guides(fill=guide_legend(title="  SD\nValue"))




heatmap.long %>%
  ggplot(aes(x = HPH, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors=color.scheme, limits = c(-2.3,2.3)) + 
  xlab("Processing Energy Level ") +
  ylab("MorFi Parameter") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 18, vjust = 2),
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_text(size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face="bold")) +
          guides(fill=guide_legend(title="  SD\nValue")) +
    facet_grid(~Variety)



heatmap.long %>%
  ggplot(aes(x = HPH, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours=colorRampPalette(c("#a50026","#f46d43", "#fee08b","white","white","#a6d96a","#66bd63","#006837"))(5), limits = c(-2.5,2.5)) + 
  xlab("Processing Energy Level ") +
  ylab("MorFi Parameter") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 18, vjust = 2),
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_text(size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face="bold")) +
  guides(fill=guide_legend(title="  SD\nValue")) +
  facet_grid(~Variety)



heatmap.long %>%
  ggplot(aes(x = HPH, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours=colorRampPalette(c("#a50026","#f46d43", "#fee08b","white","white","#a6d96a","#66bd63","#006837"))(3), limits = c(-2.5,2.5)) + 
  xlab("Processing Energy Level ") +
  ylab("MorFi Parameter") +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 18, vjust = 2),
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_text(size = 14, face="bold"),
        strip.text.x = element_text(size = 12, face="bold")) +
  guides(fill=guide_legend(title="  SD\nValue")) +
  facet_grid(~Variety)






###################### HEATMAP + Dendrogram ##########################

sub.trees <- cut(heatmap.dend, h = 0.8)
cluster3.tree <- sub.trees$lower[[3]]
cluster3.tree

Sugargraze.factor <- filter(morfi_cnf_heatmap, Variety == "Sugargraze")

Sugargraze.mtx <- as.matrix(dplyr::select(Sugargraze.factor, -HPH, -Variety))

row.names(Sugargraze.mtx) <- Sugargraze.factor$HPH

transposed.Sugargraze.mtx <- t(Sugargraze.mtx)

heatmap.2(transposed.Sugargraze.mtx,
          Rowv = cluster3.tree,  # use the dendrogram previously calculated
          Colv = NULL, # don't mess with my columns! (i.e. keep current ordering )
          dendrogram = "row",   # only draw row dendrograms
          #breaks = seq(-2, 2, length.out = 9),  # OPTIONAL: set break points for colors
          col = color.scheme,  # use previously defined colors
          trace = "none", density.info = "none",  # remove distracting elements of plot
          xlab = "Time (mins)")














###################### Prep for Morfi_Q ##########################

morfi_cnf_heatmap <- read.csv("morfi_cnf_heatmap_all.csv")

morfi_cnf_heatmap <- morfi_cnf_heatmap %>% dplyr::group_by(Variety,Section,HPH) %>% dplyr::summarise(across(.cols=everything(),mean,na.rm=TRUE)) %>% dplyr::ungroup()

morfi_cnf_heatmap2 <- morfi_cnf_heatmap[,4:24]

morfi_cnf_heatmap2 <- as_tibble(scale(morfi_cnf_heatmap2, scale=TRUE,center = TRUE))

morfi_cnf_heatmap <- bind_cols(morfi_cnf_heatmap[,1:3],morfi_cnf_heatmap2)

heatmap.long <- pivot_longer(morfi_cnf_heatmap,cols=4:24)

heatmap.long <- heatmap.long %>% dplyr::group_by(Variety,Section,HPH,name) %>% dplyr::summarise(value=mean(value,na.rm=TRUE))



 for(i in 1:nrow(heatmap.long)){heatmap.long$HPH[i] <- switch(as.character(heatmap.long$HPH[i]),"0"="L","0.5"="M","1"="H")}

heatmap.long$ID <- paste(heatmap.long$Variety,heatmap.long$Section,heatmap.long$HPH,sep="-")                                  



#################### MorFi Quality - Prep #################### 

qr <- read_csv("quality_ranking.csv")
qr$HPH.fac = factor(qr$HPH, levels=c("L", "M", "H"))
qr$ID = str_c(qr$variety,'-',qr$section,'-',qr$HPH)

Q1 <- qr %>%  filter(str_detect(quality_def,"Q1"))


morfi_df <- as.data.frame(morfi_df)

morfi_df2 <- morfi_df
morfi_df2 <- tibble::rownames_to_column(morfi_df2, "ID")
morfi_df2

morfi_Q <- merge(Q1, morfi_df2, by = "ID")
morfi_Q <- subset(morfi_Q,select = -c(variety,section,HPH,quality_def,HPH.fac))

write.csv(morfi_Q,"morfi_Q.csv", row.names = FALSE)



#################### MorFi Quality - ALL ####################                                    


morfi_quality <- join(heatmap.long,morfi_Q[,1:2],by="ID")

morfi_quality$HPH <- factor(morfi_quality$HPH, levels = c("L","M","H"))


ggplot(morfi_quality, 
       aes(x=score,y=value,color=name)) +
  geom_point(aes(alpha=abs(value)/2)) +
  geom_line(aes(alpha=abs(value)/2)) +
  scale_x_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90))


  




#################### MorFi Quality - Facet ####################


morfi_quality$HPH.fac <- factor(morfi_quality$HPH, levels=c("L", "M", "H"))


ggplot(morfi_quality, 
       aes(x=score,y=value,color=name)) +
  geom_point(aes(alpha=abs(value)/2)) +
  geom_line(aes(alpha=abs(value)/2)) +
  scale_x_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows=vars(HPH.fac)) +
  theme_bw()



ggplot(morfi_quality, 
       aes(x=value,y=score,color=name)) +
  geom_text(aes(x=value,y=score, label=name), angle = 0) +
  scale_y_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()


ggplot(morfi_quality, 
       aes(x=value,y=score,color=name)) +
  geom_text_repel(aes(x=value,y=score, label=name), angle = -45, box.padding = 0, max.overlaps = Inf, point.size = NA) +
  scale_y_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  facet_grid(cols=vars(HPH.fac))



### Less labels

ggplot(morfi_quality, 
       aes(x=value,y=score,color=name)) +
  geom_text_repel(aes(x=value,y=score, label=name), angle = -45, box.padding = 0, max.overlaps = 10, point.size = NA) +
  scale_y_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.text=element_text(size=22),
        legend.title=element_text(size=28)) +
  theme_bw() +
  facet_grid(cols=vars(HPH.fac))




### FIREWORK graph

ggplot(morfi_quality, 
       aes(x=value,y=score,color=name)) +
  # geom_point(aes(alpha=abs(value)/2)) +
  # geom_line(aes(alpha=abs(value)/2)) +
  geom_text_repel(aes(x=value,y=score, label=name), angle = 0, box.padding = 0.5, max.overlaps = Inf) +
  scale_y_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw() +
  facet_grid(cols=vars(HPH.fac))






#################### MorFi Quality - Example ####################


morfi_Q2 <- read_csv("morfi_Q2.csv")

morfi_Q2$ID <- factor(morfi_Q2$ID, levels = c("Graingrass-Leaf-L","Graingrass-Leaf-M","Graingrass-Leaf-H"))

ggplot(morfi_Q2, aes(x=score, y=value, group = parameter)) +
  geom_line(aes(colour = factor(parameter))) +
  geom_point(aes(colour = factor(parameter))) +
  xlab("Quality Score") +
  ylab("MorFi Parameter SD") +
  theme_bw() +
  labs(colour="MorFi Parameter") +
  facet_grid(rows=vars(ID,score))+
  geom_text_repel(aes(label=parameter), 
                  size = 3, 
                  fontface = "bold", 
                  hjust = -2, 
                  direction = "y", 
                  min.segment.length = 1)



# ggplot(morfi_quality, 
#        aes(x=score,y=value,color=name)) +
#   geom_point(aes(alpha=abs(value)/2)) +
#   geom_line(aes(alpha=abs(value)/2)) +
#   scale_x_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_grid(HPH~., margins = TRUE)
# 
# 
# ggplot(morfi_quality, 
#        aes(x=score,y=value,color=name)) +
#   geom_point(aes(alpha=abs(value)/2)) +
#   geom_line(aes(alpha=abs(value)/2)) +
#   scale_x_continuous(breaks=morfi_quality$score,labels=morfi_quality$ID) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_grid(rows=vars(name))





