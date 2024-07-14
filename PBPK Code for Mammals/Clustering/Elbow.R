rm(list = ls())

df <- read.csv('C:/XXX/OneDrive - Syngenta/AI/Propiconazole/PBK_Modeling/in vitro/4clustering_1007.csv')
kmax <- 5

row.names(df) <- df$Compound
df[,1]        <- NULL
dd            <- dist(scale(df), method = "euclidean")
hc  <- hclust(dd, method = "ward.D2")

# Compute the total within-cluster sum of square (WSS) for different number of clusters
wss <- sapply(1:kmax, function(k) {
  hc.k <- cutree(hc, k)
  sum(sapply(1:k, function(i) {
    cluster_points <- df[hc.k == i, ]
    sum(dist(cluster_points)^2)
  }))
})

# Plot the total within-cluster sum of square (WSS) vs. the number of clusters
plot(1:kmax, wss, type = "b", xlab = "Number of clusters", ylab = "Within groups sum of squares")

