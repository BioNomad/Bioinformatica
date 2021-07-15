Clustering
================

Clustering, like the name implies, is a way to group data points to find
patterns. Clustering is type of unsupervised learning - where patterns
are discovered without labeled data
([Wikipedia](https://en.wikipedia.org/wiki/Cluster_analysis)). Here we
are going to talk about k-means clustering and hierarchical clustering.

## K-means Clustering

In k-means clustering, you designate how many clusters are to be found
in the data ([Wikipedia](https://en.wikipedia.org/wiki/K-)). The k-means
algorithm will randomly assign this many “cluster centers” and the rest
of the data points will be assigned to those clusters by minimizing
Euclidean distance ([Wikipedia](https://en.wikipedia.org/wiki/K-)). This
process goes through multiple iterations until the variance within those
clusters is minimized ([Wikipedia](https://en.wikipedia.org/wiki/K-)).
Let’s try it in R\!

``` r
load("./lgg.rda")
library(limma)
library(factoextra)

#now let's normalize!
log.trans <- log2(lgg$ExpressionData + 1)
norm.data <- normalizeQuantiles(log.trans)
#let's just look at the top 500 genes with the highest variance
vars <- apply(norm.data,1,var)
select = names(vars[order(vars,decreasing = T)][1:500])
top.vars = norm.data[select,]

#how many clusters should we call?
fviz_nbclust(top.vars,kmeans) 
```

![](clustering_files/figure-gfm/kmean-1.svg)<!-- -->

``` r
#alright let's get those clusters and visualize them
k_res <- kmeans(top.vars,2,nstart = 25)
fviz_cluster(k_res, geom = "point", data = top.vars) + ggtitle("k = 2")
```

![](clustering_files/figure-gfm/kmean-2.svg)<!-- -->

So as we saw, the recommended number of clusters for these data was 2.
There also doesn’t seem to be a great separation of clusters. This
highlights a major caveat of k-means clustering - it is highly dependent
on the shape of your data ([towards data
science](https://towardsdatascience.com/k-means-clustering-algorithm-applications-evaluation-methods-and-drawbacks-aa03e644b48a)).

## Hierarchal Clustering

So what about hierarchical clustering? Hierarchical clustering works by
building well a hierarchy of clusters. They can either be generated via
an agglomerative method or divisive
([Wikipedia](https://en.wikipedia.org/wiki/Hierarchical_clustering)).
Agglomerative clustering generates clusters from the bottom up, where
each data point starts out as its own cluster
([Wikipedia](https://en.wikipedia.org/wiki/Hierarchical_clustering)).
Divisive clustering takes all data points as one cluster and then splits
them into multiple clusters
([Wikipedia](https://en.wikipedia.org/wiki/Hierarchical_clustering)).
Agglomerative clusters are created by minimizing
([Wikipedia](https://en.wikipedia.org/wiki/Hierarchical_clustering)):

  - The maximum distance between elements of each cluster (also called
    complete-linkage clustering)

  - The minimum distance between elements of each cluster (also called
    single-linkage clustering):

  - The mean distance between elements of each cluster (also called
    average linkage clustering, used e.g. in UPGMA):

  - The sum of all intra-cluster variance.

  - The increase in variance for the cluster being merged (Ward’s
    method)

  - The probability that candidate clusters spawn from the same
    distribution function (V-linkage).

Let’s try Ward’s clustering in R\!

``` r
#let's create a distance matrix for
distances <- get_dist(t(top.vars),method = "euclidean")
#run ward's clustering
hc <- hclust(distances,method = "ward.D2")
#plot the dendrogram
plot(hc,labels = FALSE)
abline(h=450)
```

![](clustering_files/figure-gfm/wards-1.svg)<!-- -->

``` r
#let's cut the dendrogram to create 3 clusters
clusters <- cutree(hc,h=450)
#now let's plot our clusters!
fviz_cluster(list(data = t(top.vars), cluster = clusters),geom = "point")+
  ggtitle("Ward's Clustering n=3")
```

![](clustering_files/figure-gfm/wards-2.svg)<!-- -->

Here we can see that this time we need a distance matrix. We just need
to be sure to make sure that when we calculate those distances using the
`get_dist()` function, that the points we want to compare are the rows
of the data frame. So here we had to use the transpose function, `t()`,
to make sure patients are the rows, not the columns. We also generate a
dendrogram - a visualization of the hierarchy. We then “cut” this
dendrogram at a given height to create the clusters.

## References

1.  <https://en.wikipedia.org/wiki/Cluster_analysis>

2.  <https://en.wikipedia.org/wiki/K-means_clustering>

3.  <https://towardsdatascience.com/k-means-clustering-algorithm-applications-evaluation-methods-and-drawbacks-aa03e644b48a>

4.  <https://en.wikipedia.org/wiki/Hierarchical_clustering>
