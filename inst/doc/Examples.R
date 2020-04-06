## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- collapse=TRUE-----------------------------------------------------------
# Package `FunChisq' must have been installed.

plot.patterns <- function(x, y, z, res)
{
  k = length(unique(res$clabels))
  mar = c(2.5,3,4,1)
  par(mar=mar, mgp=c(3,1,0)-c(1.5,0.5,0), lwd=2)
  
  col <- "limegreen"
  labelcol <- colorRampPalette(c("black", col))
  plot(x, y, main="Original data", col=labelcol(k)[res$clabels], 
       pch=19, cex.axis=0.8, cex=0.7)
  abline(v=res$grid[[1]], h=res$grid[[2]], col="black", lty="dotted")
  tab <- as.matrix(table(-res$D[, 2], res$D[, 1]))
  par(xpd=TRUE)
  legend("top", legend = c(min(res$clabels):max(res$clabels)), pch=19,
         col = labelcol(k), bty = "n", horiz = TRUE, inset = c(0,-0.2))
  par(xpd=FALSE)
  FunChisq::plot_table(
    tab, xlab="x discretized", ylab="y discretized",
    col=col, main="Discretized data", highlight="none", mar = mar)
  
  col <- "brown3"
  labelcol <- colorRampPalette(c("black", col))
  plot(y, z, main="Original data", col=labelcol(k)[res$clabels], 
       pch=19, cex.axis=0.8, cex=0.7)
  abline(v=res$grid[[2]], h=res$grid[[3]], col="black", lty="dotted")
  tab <- as.matrix(table(-res$D[, 3], res$D[, 2]))
  par(xpd=TRUE)
  legend("top", legend = c(min(res$clabels):max(res$clabels)), pch=19, 
         col = labelcol(k), bty = "n", horiz = TRUE, inset = c(0,-0.2))
  par(xpd=FALSE)
  FunChisq::plot_table(
    tab, xlab="y discretized", ylab="z discretized",
    col=col, main="Discretized data", highlight="none", mar = mar) 
  
  col <- "dodgerblue"
  labelcol <- colorRampPalette(c("black", col))
  plot(x, z, main="Original data", col=labelcol(k)[res$clabels],
       pch=19, cex.axis=0.8, cex=0.7)
  abline(v=res$grid[[1]], h=res$grid[[3]], col="black", lty="dotted")
  tab <- as.matrix(table(-res$D[, 3], res$D[, 1]))
  par(xpd=TRUE)
  legend("top", legend = c(min(res$clabels):max(res$clabels)), pch=19,
         col = labelcol(k), bty = "n", horiz = TRUE, inset = c(0,-0.2))
  par(xpd=FALSE)
  FunChisq::plot_table(
    tab, xlab="x discretized", ylab="z discretized",
    col=col, main="Discretized data", highlight="none", mar = mar) 
}

## ---- out.width="40%", fig.show="hold", fig.cap="Example 1. Nonlinear curves using k-means clustering with a fixed number of clusters."----
require(GridOnClusters)
x = rnorm(50)
y = sin(x)
z = cos(x)
data = cbind(x, y, z)
res = discretize.jointly(data, k=3) # using a specified k
plot.patterns(x, y, z, res)

## ---- out.width="40%", fig.show="hold", fig.cap="Example 2. Using a range for the number of k-means clusters"----
 x = rnorm(100)
 y = log1p(abs(x))
 z = ifelse(x >= -0.5 & x <= 0.5, 0, 1) + rnorm(100, 0, 0.1)
 data = cbind(x, y, z)
 res = discretize.jointly(data, k=c(2:3)) # using a range of k
 plot.patterns(x, y, z, res)

## ---- out.width="40%", fig.show="hold", fig.cap="Example 3. Using the partition around medoids clustering method."----
 # using a clustering method other than k-means
 x = rnorm(100)
 y = log1p(abs(x))
 z = sin(x)
 data = cbind(x, y, z)

 # pre-cluster the data using partition around medoids (PAM)
 cluster_label = cluster::pam(x=data, diss = FALSE, metric = "euclidean", k = 4)$clustering
 res = discretize.jointly(data, cluster_label = cluster_label)
 plot.patterns(x, y, z, res)

