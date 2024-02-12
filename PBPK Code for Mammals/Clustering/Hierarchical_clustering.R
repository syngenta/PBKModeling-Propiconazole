library(factoextra)
library(dendextend)
library(ggsignif)
library(stats)
# reference : 
# https://www.datanovia.com/en/lessons/examples-of-dendrograms-visualization/
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/
# http://agroninfotech.blogspot.com/2020/06/visualizing-clusters-in-r-hierarchical.html
rm(list = ls())


# Compute distances and hierarchical clustering
df <- read.csv('C:/XXX/in vitro/4clustering_0811.csv')
df <- read.csv('C:/XXX/in vitro/4clustering_1007.csv')
#df <- read.csv('C:/XXX/in vitro/4clustering_0102.csv')
row.names(df) <- df$Compound
df[,1]        <- NULL
dd            <- dist(scale(df), method = "euclidean")

hc <- hclust(dd, method = "ward.D2")
fviz_dend(hc, cex = 0.5)
fviz_dend(hc, cex = 0.5, horiz = TRUE)

fviz_dend(hc, k = 7,                 # Cut in four groups
          cex = 0.5,                 # label size
          horiz = TRUE,
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800","#FC4E07", 'black', 'blue','orange'),
          color_labels_by_k = TRUE,  # color labels by groups
          ggtheme = theme_gray()     # Change theme
)+ theme_bw()

tiff("C:/XXX/clustering_1008.tiff", units="in", width=6, height=5, res=600,compression = "lzw")

