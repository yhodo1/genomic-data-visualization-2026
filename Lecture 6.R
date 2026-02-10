data <- read.csv("~/Desktop/Visium-IRI-ShamR_matrix.csv.gz")
dim(data)
pos <- data[,c('x', 'y')]
  
rownames(gexp) <- data [,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

#normalize
totgexp = rowSums(gexp)
mat <- log10(gexp/totgexp * 1e6 + 1)

#PCA
pcs <- prcomp(mat)
df <- data.frame(pcs$x, pos)
ggplot(df, aes(x=x, y=y, col=PC1)) + geom_point(cex=0.5)

#kmeans
clusters <- as.factor(kmeans(pcs$x[,2:2], centers=2)$cluster)
df <- data.fram(pcs$x, pos, clusters)
ggplot(df, aes(x=PC1, y=PC2, col=clusters)) + geom_point(cex=0.5)
ggplot(df, aes(x=x, y=y, col=clusters)) + gemon_point(cex=0.5)

#differential expression
?t.test
? wilcox.test
'Nphs2' %in% colnames(mat)
df <- data.frame(pcs$x, pos, clusters, gene=mat[,'Nphs2'])
ggplot(df, aes(x=PCS, y=PC6, col=gene)) + geom_point(cex=0.5, alpha=0.5)

#between cells of cluster 2 vs all others 
clusterofintersest <- names(clusters)[clusters ==2]
othercells <- names(clusters)[clusters !=2]
x1 <- mat[othercells, 'Nphs2']
wilcox.test(x1, x2, alternative = ' two.sided') # different expressed
wilcox.test(x1, x2, alternative= 'greater') #differentially upregulated in cluster 
wilcox.test(x2, x1, alternative= 'greater') # not differentially upregulated in other cells
wilcox.test(x2, x1, alternative= 'greater')$p.value #grab p-value

#loop
?apply
?sapply
?lapply
out <- sapply(colnames(mat), function(gene){
  x1 <- mat[clusterofinterest, gene]
  x2 <- mat[othercells, gene]
  wilcox.test(x1, x2, alternative= 'greater')$p.value
})
head(sort(out))
df <- data.frame(pcs$x, pos, clusters, gene=mat[, 'Cd24a'])
ggplot(df, aes(xPC1, y=PC2, col=gene)) + geom_point(cex=0.5, alpha=0.5)
