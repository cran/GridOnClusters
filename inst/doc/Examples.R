## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- out.width="40%", fig.show="hold", fig.cap="Example 1. Nonlinear curves using $k$-means clustering with a range for the number of clusters."----
require(GridOnClusters)
require(FunChisq)

plot.patterns <- function(x, y, z, res)
{
 par(mar=c(2.5,3,3,1), mgp=c(3,1,0)-c(1.5,0.5,0))
 
 col <- "green3"
 plot(x, y, main="Original data", col=col)
 abline(v=res$grid[[1]], h=res$grid[[2]], col=col, lty="dotted")
 tab <- as.matrix(table(-res$D[, 2], res$D[, 1]))
 plot_table(tab, xlab="x discretized", ylab="y discretized",
            col=col, main="Discretized data", highlight="none")
 
 col <- "red3"
 plot(y, z, main="Original data", col=col)
 abline(v=res$grid[[2]], h=res$grid[[3]], col=col)
 tab <- as.matrix(table(-res$D[, 3], res$D[, 2]))
 plot_table(tab, xlab="y discretized", ylab="z discretized",
            col=col, main="Discretized data", highlight="none") 
 
 col <- "royalblue"
 plot(x, z, main="Original data", col=col)
 abline(v=res$grid[[1]], h=res$grid[[3]], col=col)
 tab <- as.matrix(table(-res$D[, 3], res$D[, 1]))
 plot_table(tab, xlab="x discretized", ylab="z discretized",
            col=col, main="Discretized data", highlight="none") 
}
   
# using a specified \code{k}
x = rnorm(40)
y = sin(x)
z = cos(x)
data = cbind(x, y, z)
res = discretize.jointly(data, k=5)
plot.patterns(x, y, z, res)

## ---- out.width="40%", fig.show="hold", fig.cap="Example 2. Using fixed number of $k$-means clusters"----
 # using a range of 'k'
 n = 40
 x = rnorm(n)
 y = log1p(abs(x))
 z = ifelse(x >= -.7 & x <= .7, 0, 1) + rnorm(n, 0, 0.175)
 data = cbind(x, y, z)
 res = discretize.jointly(data, k=c(2:3))
 plot.patterns(x, y, z, res)

## ---- out.width="40%", fig.show="hold", fig.cap="Example 3. Using the partition around medoids clustering method."----
 # using an alternate different clustering scheme
 library(cluster)
 n = 40
 x = rnorm(n)
 y = log1p(abs(x))
 z = sin(x)
 data = cbind(x, y, z)

 # pre-cluster the data using partition around medoids (PAM)
 cluster_label = pam(x=data, diss = FALSE, metric = "euclidean", k = 5)$clustering
 res = discretize.jointly(data, cluster_label = cluster_label)
 plot.patterns(x, y, z, res)
