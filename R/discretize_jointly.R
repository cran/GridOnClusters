# discretize.jointly.R
#
# Created by Sajal Kumar
# Modified by Jiandong Wang
# Copyright (c) NMSU Song lab

#' @import Rcpp
#' @useDynLib GridOnClusters
#' 
#' @importFrom cluster silhouette
#' @importFrom stats kmeans
#' @importFrom fossil adj.rand.index
#' @importFrom stats dist
#' @importFrom dqrng dqsample
#' @importFrom Rdpack reprompt
#' @importFrom mclust Mclust mclustBIC
#' @importFrom Ckmeans.1d.dp MultiChannel.WUC ahist
#  @importFrom Telescope telescope


#' @title Discretize Multivariate Continuous Data by Cluster-Preserving Grid
#'
#' @description
#'  Discretize multivariate continuous data using a grid that captures 
#' the joint distribution via preserving clusters in original data
#'
#' @param data a numeric matrix for multivariate data or a numeric vector for 
#'   univariate data. In case of a matrix, columns are continuous variables; 
#'   rows are observations.
#'
#' @param k either an integer, a integer vector, 
#'   or \code{Inf}, specifying the number of clusters. 
#'   The default is a vector of integers from 2 to 10. 
#'   If \code{k} is a single number, \code{data} will 
#'   be grouped into into exactly \code{k} clusters. 
#'   If \code{k} is an integer vector, an optimal 
#'   \code{k} is chosen among the integers. If \code{k} 
#'   is set to \code{Inf}, an optimal \code{k} is 
#'   chosen from 2 to \code{nrow(data)}. If 
#'   \code{cluster_label} is specified, 
#'   \code{k} is ignored.
#'
#' @param min_level an integer or an integer vector, to specify the minimum number of levels 
#' along each dimension. If a vector of size \code{ncol(data)}, then each element
#' will be mapped 1:1 to each dimension in order. If an integer, then all dimensions
#' will have the same minimum number of levels.
#'
#' @param max_level an integer or an integer vector, to specify the maximum
#'   number of levels along each dimension. It works in the 
#'   same way as \code{min_level}. \code{max_level} will
#'   be set to the smaller between number of compressed zones and itself,
#'   if \code{grid_method} is a likelihood approach or 
#'   "DP Compressed majority".
#'
#' @param cluster_method a character string to specify a clustering
#'  method to be used. Ignored if \code{cluster_label} is not \code{NULL}.
#'  We offer three build-in options:
#'  
#' \code{"Ball+BIC"} (default) uses \code{mclust::Mclust} 
#'   (\code{modelNames = "VII"} for 2-D or higher dimensions; 
#'   \code{"V"} for 1-D) to cluster \code{data} and 
#'   BIC score to select number of clusters. 
#' 
#' \code{"kmeans+silhouette"} uses k-means to cluster \code{data} and the average 
#' Silhouette width to select number of clusters.
#' 
#' \code{"PAM"} uses the algorithm partition around medoids to perform clustering.
#' 
#' @param grid_method a character string to specify a grid 
#' discretization method. Default: 
#' \code{"DP approx likelihood 1-way"}. The methods 
#' can be roughly separate into three different categories: 
#' by cluster likelihood, by density, and by SSE (Sum of Squared Errors).
#' See Details for more information.
#'  
#' @param eval_method a character string to 
#'   specify a method to evaluate quality 
#'   of discretized data.
#'  
#' @param cluster_label a vector of labels for each data point or 
#'   observation. It can be class labels on the input \code{data} for 
#'   supervised learning; it can also be cluster labels for 
#'   unsupervised learning. If \code{NULL} (default), clustering 
#'   is performed to obtain labels.
#'
#' @param cutoff a numeric value. A grid line is added only when the 
#'   quality of the line is not smaller than \code{cutoff}. 
#'   It is applicable only to \code{grid_method} \code{"DP"} or 
#'   \code{"DP Compressed majority"}.
#' 
#' @param entropy a logical to chose either entropy 
#'   (\code{TRUE}) or likelihood (\code{FALSE}, default).
#'   
#' @param noise a logical to apply jitter noise to original 
#'   data if \code{TRUE}. Default: \code{FALSE}. 
#'   It is only applicable 
#'   to \code{cluster_method} \code{"Ball+BIC"}. 
#'   When data contain many duplicated values, 
#'   adding noise can help \code{Mclust} clustering.  
#'  
#' @param dim_reduction a logical to turn on/off 
#'   dimension reduction. Default: \code{FALSE}.
#'   
#' @param scale a logical to specify linear 
#'   scaling of the variable in each dimension
#'   if \code{TRUE}. Default: \code{FALSE}.
#'   
#' @param variance a numeric value to specify 
#'   noise variance to be added to the data

#' @param nthread an integer to specify number 
#'   of CPU threads to use. Automatically adjusted
#'   if invalid or exceeding available cores.
#'
#'
#' @note
#' The default \code{grid_method} is changed
#' from \code{"Sort+Split"} \insertCite{Jwang2020BCB}{GridOnClusters} (up to released package version 0.1.0.2) 
#' to \code{"DP approx likelihood 1-way"} (since version 0.3.2), 
#' representing a major improvement.
#' 
#' 
#' @details 
#' 
#' The function implements both published algorithms described in 
#' \insertCite{Jwang2020BCB}{GridOnClusters} and new algorithms for
#' multivariate discretization.
#' 
#' The included grid discretization methods can be summarized into three categories: 
#' 
#' * By Density
#'
#'     - \code{"Sort+split"} \insertCite{Jwang2020BCB}{GridOnClusters} 
#'     sorts clusters by mean in each dimension. It then
#'  splits consecutive pairs only if the sum of error rate of each cluster is
#'  less than or equal to 50%. It is possible that no grid line will be added
#'  in a certain dimension. The maximum number of lines is the number of 
#'  clusters minus one.
#'  
#' * By SSE (Sum of Squared Errors)
#'     - \code{"MultiChannel.WUC"} splits each dimension by weighted with-in cluster
#'  sum of squared distances by \code{Ckmeans.1d.dp::MultiChannel.WUC()}. Applied in 
#'  each projection on each dimension. The channel of each point is defined by 
#'  its multivariate cluster label.
#'  
#'     - \code{"DP"} orders labels by data in each dimension and then cuts data
#'  into a maximum of \code{max_level} bins. It evaluates the quality of each 
#'  cut to find a best number of bins.
#'  
#'     - \code{"DP Compressed majority"} orders labels by data in each dimension. 
#'  It then compresses labels neighbored by the same label to avoid 
#'  discretization within consecutive points of the same cluster label, so as to 
#'  greatly reduce runtime of dynamic programming. Then it cuts data into 
#'  a maximum of \code{max_level} bins, and it evaluates the quality of 
#'  each cut by the majority of data to find a best number of bins.
#'  
#' * By cluster likelihood
#'     - \code{"DP exact likelihood"} orders labels by data in each dimension.
#'  It then compresses labels neighbored by the same label to avoid 
#'  discretization within consecutive points of the same cluster label, 
#'  so as to greatly reduce runtime of dynamic programming. 
#'  Then cut the data into a maximum of \code{max_level} bins.
#'  
#'     - \code{"DP approx likelihood 1-way"} is a sped-up version of the 
#'  \code{"DP exact likelihood"} method, but it is not always optimal.
#' 
#'     - \code{"DP approx likelihood 2-way"} is a bidirectional variant of the 
#'  \code{"DP approx likelihood"} method. It performs approximate dynamic 
#'  programming in both the forward and backward directions and selects 
#'  the better of the two results. This approach provides additional robustness 
#'  compared to the one-directional version, but optimality is not always achieved.
#'  
#' 
#' @return
#'
#' A list that contains four items:
#' \item{\code{D}}{a matrix of discretized values from original \code{data}.
#' Discretized values are one(1)-based.}
#'
#' \item{\code{grid}}{a list of numeric vectors of decision boundaries for each variable/dimension.}
#'
#' \item{\code{clabels}}{a vector of cluster labels for each observation in \code{data}.}
#'
#' \item{\code{csimilarity}}{a similarity score between clusters from joint discretization
#' \code{D} and cluster labels \code{clabels}. The score is the adjusted Rand index.}
#' 
#' @examples
#' # using a specified k
#' x = rnorm(100)
#' y = sin(x)
#' z = cos(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=5)$D
#'
#' # using a range of k
#' x = rnorm(100)
#' y = log1p(abs(x))
#' z = tan(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=c(3:10))$D
#' 
#' # using k = Inf
#' x = c()
#' y = c()
#' mns = seq(0,1200,100)
#' for(i in 1:12){
#'   x = c(x,runif(n=20, min=mns[i], max=mns[i]+20))
#'   y = c(y,runif(n=20, min=mns[i], max=mns[i]+20))
#' }
#' data = cbind(x, y)
#' discretized_data = discretize.jointly(data, k=Inf)$D
#'
#' # using an alternate clustering method to k-means
#' library(cluster)
#' x = rnorm(100)
#' y = log1p(abs(x))
#' z = sin(x)
#' data = cbind(x, y, z)
#'
#' # pre-cluster the data using partition around medoids (PAM)
#' cluster_label = pam(x=data, diss = FALSE, metric = "euclidean", k = 5)$clustering
#' discretized_data = discretize.jointly(data, cluster_label = cluster_label)$D
#'
#' @references
#' \insertAllCited{}
#' 
#' @author 
#' 
#' Jiandong Wang, Sajal Kumar, and Mingzhou Song
#' 
#' @seealso
#'
#' See \link[Ckmeans.1d.dp]{Ckmeans.1d.dp} for discretizing univariate continuous data.
#'
#' @md
#'
#' @export
discretize.jointly = function(
    data, k=c(2:10), min_level = 1, max_level= 100,
    cluster_method = c(
      "Ball+BIC", 
      "kmeans+silhouette", 
      "PAM" # , 
      # "telescope",
      # "Ball+Tele+BIC"
    ),
    grid_method = c(
      "DP approx likelihood 1-way", 
      "DP approx likelihood 2-way",
      "DP exact likelihood", 
      "DP Compressed majority", 
      "DP", 
      "Sort+split", 
      "MultiChannel.WUC"),
    eval_method = c(
      "ARI", "purity", "upsllion", "CAIR"
    ),
    cluster_label=NULL, cutoff=0.00, 
    entropy=FALSE,
    noise = FALSE, 
    dim_reduction = FALSE, 
    scale=FALSE,
    variance = 0.5, 
    nthread=1
)
{

  # check if data provided is a matrix
  # if( !("matrix" %in% class(data)) && !("data.frame" %in% class(data))){
  #  stop("'data' must be a matrix or data.frame.")
  # }
  
  if(inherits(data, c("numeric", "integer"))) {
    # Check if a numeric or integer vector:
    data <- matrix(data, ncol = 1)
  } else if(! inherits(data, c("matrix", "data.frame"))) {
    # Check if data provided is a matrix or data frame
    stop("'data' must be a matrix or data.frame.")
  }
  
  # check if all columns in 'data' are numeric or integer
  dim_class = apply(data, 2, class)
  if(!all(dim_class == "numeric" | dim_class == "integer")){
    stop("All columns in 'data' should be numeric or integer.")
  }

  # 'data' should have at least 10 observations
  dim_data = dim(data)
  #if(dim_data[1] < 10){
  #  stop("'data' should have at least 10 observations,")
  #}

  # 'cluster_label' should either be null or integers (or numeric) matching nrow(data)
  if(!is.null(cluster_label) && length(cluster_label) != nrow(data)){
    stop("'cluster_label' should either be null or a vector with nrow(data) elements.")
  }

  if(!is.null(cluster_label) && !class(cluster_label) %in% c("numeric","integer")){
    stop("'cluster_label' should be either null or a numeric/integer vector.")
  }
  
  # if "min_level" is a list or a integer
  if(length(min_level)==1){
    min_level = rep(min_level, ncol(data))
  } else if(length(min_level) > 1 && length(min_level) != ncol(data)){
    stop("'min_level' should either be an integer or a vector of size ncol(data)")
  }
  
  # if "max_level" is a list or a integer
  if(length(max_level)==1){
     max_level = rep(max_level, ncol(data))
  } else if(length(max_level) > 1 && length(max_level) != ncol(data)){
     stop("'min_level' should either be an integer or a vector of size ncol(data)")
  }
    
  # switching order, checking min_level < k after having a k
  
  # 'min_level' should smaller then the max of k
  if(max(k) < max(min_level)){
    k = c(min_level,min_level+5)
    warning("'min_level' should be in the range of k, k has been adapted")
   } 
  
  # the maximum of k range should be less then the number of point
  if(all(is.finite(k))){
    k = k[which(k<nrow(data))]  
  }
  
  #scale data
  data_scaled = data
  if(scale == TRUE){
     data_scaled = scale(data)
  }
  #check method
  cluster_method = match.arg(cluster_method)
  grid_method = match.arg(grid_method)
  eval_method = match.arg(eval_method)
  
  # if no cluster labels and methods are supplied, default to mclust (VII or V for 1_D)
  if(is.null(cluster_label)){

    cluster_info =  cluster(data_scaled, k, cluster_method, noise)

  } else { # cluster labels are supplied

    cluster_label = as.numeric(as.factor(cluster_label)) # making labels consecutive

    # # find medoids
    # centers = data[medoids(D = dist(data), cl = cluster_label),,drop=FALSE]

    cluster_info = list(#centers = as.matrix(centers),
                        clusters = cluster_label-1,
                        data = as.matrix(data), cluster_method = "user supply")
  }

  #browser()
  if(grid_method == "MultiChannel.WUC"){  
    # prepare the weights
    # kstar = length(unique(cluster_info$clusters)) # MS 2022-01-17
    clusters.factor <- factor(cluster_info$clusters) # MS 2022-01-17
    kstar <- nlevels(clusters.factor)                # MS 2022-01-17
    weight = matrix(1/nrow(data), nrow = nrow(data), ncol = kstar)
    for (i in c(1:nrow(data))) {
      # weight[i, cluster_info$clusters[i]+1] = nrow(data)/(nrow(data)+1) # MS 2022-01-17
       weight[i, clusters.factor[i]] = nrow(data)/(nrow(data)+1) # MS 2022-01-17
    }
    grid_lines = vector("list", length = ncol(data))
    # for n dim
    for(dim in c(1:ncol(data))){
      result = Ckmeans.1d.dp::MultiChannel.WUC(x = data[,dim], y = weight, k = c(min_level[dim]:kstar))
      result$size = as.vector(table(result$cluster))
      res <- structure(result,class = "Ckmeans.1d.dp")
      
      res.hist = Ckmeans.1d.dp::ahist(res, data = data[,dim], style = "midpoints", plot = FALSE)
      
      lines = unique(res.hist$breaks)
      grid_lines[[dim]] = lines[-c(1,length(lines))]
    }
  }

  else if(grid_method == "Sort+split"){
    # get grid lines
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), nrow(data), 
                     ncol(data), min_level, max_level, "Density", cutoff, entropy, 
                     FALSE, nthread)
    # filter grid lines
    #grid_lines = lapply(grid_lines, function(i){ 
    #  mx = max(i) 
    #  return(i[-which(i == mx)])
    #  })
  }
  else if(grid_method == "DP"){
     # get grid lines
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), 
                     nrow(data), ncol(data), min_level, max_level, "DP majority",
                     cutoff, entropy, FALSE, nthread)
     
  }
  else if(grid_method == "DP Compressed majority"){
     # get grid lines
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), 
                     nrow(data), ncol(data), min_level, max_level, "DP compressed majority",
                     cutoff, entropy, FALSE, nthread)
     
  }
  else if(grid_method == "DP exact likelihood"){
     # get grid lines
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), 
                     nrow(data), ncol(data), min_level, max_level, "DP exact likelihood", 
                     0, entropy, dim_reduction, nthread)
  }
  else if(grid_method == "DP approx likelihood 1-way"){
     # get grid lines
     #browser()
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), 
                     nrow(data), ncol(data), min_level, max_level, "DP approx likelihood",
                     0, entropy, dim_reduction, nthread)
     #return(grid)
  }
  else if(grid_method == "DP approx likelihood 2-way"){
     # get grid lines
     #browser()
     grid = findgrid(cluster_info, length(unique(cluster_info$clusters)), 
                     nrow(data), ncol(data), min_level, max_level, "DP approx likelihood two way",
                     0, entropy, dim_reduction, nthread)
     #return(grid)
  }
  names(grid) <- c("lines", "table", "discr", "BIC", "purity", "lines_reduc", 
                   "table_reduc", "discr_reduc", "dimensional_bic", "BIC_reduc",
                   "purity_reduc", "line_removed", "upsilon", 
                   "distance_to_median", "distance_to_mean")

  if(FALSE){
  # discretize data
  discr_data = discretize_data(data, grid$lines)

  if(eval_method == "ARI"){
     # compute adjusted random index
     ndim_cluster_dist = discr_data[,1]
     if(ncol(discr_data) > 1){
        for(i in 2:ncol(discr_data)){
           ndim_cluster_dist = paste0(ndim_cluster_dist,",",discr_data[,i])
        }
     }
     #browser()
     block_ID = as.factor(ndim_cluster_dist)
     if(nlevels(block_ID) > 1){
        cluster_similarity = mclust::adjustedRandIndex(as.numeric(block_ID)-1, cluster_info$clusters)
     }
     else{
        cluster_similarity = 0
     }
  }
  else if(eval_method == "purity"){
     # compute adjusted random index
     ndim_cluster_dist = discr_data[,1]
     if(ncol(discr_data) > 1){
        for(i in 2:ncol(discr_data)){
           ndim_cluster_dist = paste0(ndim_cluster_dist,",",discr_data[,i])
        }
     }
     #browser()
     block_ID = as.factor(ndim_cluster_dist)
     #cluster_similarity = purity(cluster_info$clusters, as.numeric(block_ID)-1)
     #browser()
     cluster_similarity = ClusterPurity(as.numeric(block_ID)-1, cluster_info$clusters)
  }
  }
  
  # compute adjusted random index
  N = dim(data)[1]
  ndim = dim(data)[2]
  discr_data = grid$discr
  discr_data_after = grid$discr_reduc
  
  discr_data = matrix(unlist(discr_data), nrow=N, byrow = TRUE)
  if(dim_reduction){
     discr_data_after = matrix(unlist(discr_data_after), nrow=N, byrow = TRUE)
  }
  
  # browser()
  # compute adjusted random index for original data
  ndim_cluster_dist = discr_data[,1]
  if(ncol(discr_data) > 1){
     for(i in 2:ncol(discr_data)){
        ndim_cluster_dist = paste0(ndim_cluster_dist,",",discr_data[,i])
     }
  }
  block_ID = as.factor(ndim_cluster_dist)
  if(nlevels(block_ID) > 1){
     #cluster_similarity_o = mclust::adjustedRandIndex(as.numeric(block_ID)-1, cluster_info$clusters)
     cluster_similarity_o =  mclust::adjustedRandIndex(as.numeric(block_ID)-1, cluster_info$clusters)
  }
  else{
     cluster_similarity_o = 0
  }
  
  if(dim_reduction){
     # compute adjusted random index for removed data
     ndim_cluster_dist_r = discr_data_after[,1]
     if(ncol(discr_data_after) > 1){
        for(i in 2:ncol(discr_data_after)){
           ndim_cluster_dist_r = paste0(ndim_cluster_dist_r,",",discr_data_after[,i])
        }
     }
     block_ID = as.factor(ndim_cluster_dist_r)
     if(nlevels(block_ID) > 1){
        #cluster_similarity_r = mclust::adjustedRandIndex(as.numeric(block_ID)-1, cluster_info$clusters)
        cluster_similarity_r =  mclust::adjustedRandIndex(as.numeric(block_ID)-1, cluster_info$clusters)
     }
     else{
        cluster_similarity_r = 0
     }
  }
  else{
     cluster_similarity_r = 0
  }
  
  
  contig_table_o = matrix(unlist(grid$table), nrow=length(grid$table), byrow = TRUE)
  if(dim_reduction){
     contig_table_r = matrix(unlist(grid$table_reduc), nrow=length(grid$table_reduc), byrow = TRUE)
  }
  else{
     contig_table_r = 0
  }
  
  # calculate upsilon
  max_upsilon = (prod(apply(discr_data, 2, function(col) length(unique(col))))*(dim(data)[1]))/(2^(dim(data)[2]))
  u.res = hashed_upsilon(discr_data, TRUE)
  uval = u.res$statistic
  #urelt = uval / max_upsilon * (dim(data)[1])
  urelt = u.res$stat_relt
  pval = u.res$p.value
  #prelt = pchisq(urelt, df=1, lower.tail = FALSE, log.p = TRUE)
  prelt = u.res$p.value_relt
  
  # calculate p-val for upsilon
  #pval = 0
  #pval = pchisq(as.numeric(grid$upsilon), df = (dim(contig_table_o)[1]-1)*(dim(contig_table_o)[2]-1), 
  #              lower.tail = FALSE, log.p = TRUE)
  
  # print(cluster_similarity_o)
  # print(cluster_similarity_r)
  
  class = "GridOnClusters"
  # browser()
  res = structure(
     list(D=discr_data, D_after = discr_data_after, 
          grid=grid$lines, grid_r = grid$lines_reduc, 
          contig_table = contig_table_o, 
          contig_table_r = contig_table_r, 
          clabels=cluster_info$clusters+1, 
          dlabels=as.numeric(block_ID)-1, 
          dimensional_bic = grid$dimensional_bic, 
          csimilarity=cluster_similarity_o, 
          csimilarity_r = cluster_similarity_r, 
          BIC = grid$BIC, BIC_r = grid$BIC_reduc,
          purity = grid$purity, 
          purity_r = grid$purity_reduc, pval = pval, 
          prelt = prelt,
          upsilon = uval, upsilon_relt = urelt, 
          cluster_method = cluster_method, 
          upsilon_level = u.res$upsilon_level, 
          tsp = u.res$tsp, tsm = u.res$tsm,
          distance_to_median = grid$distance_to_median, 
          mean_dist_to_median = mean(unlist(grid$distance_to_median)),
          distance_to_mean = grid$distance_to_mean, 
          mean_dist_to_mean =  mean(unlist(grid$distance_to_mean)),
          grid_method = grid_method, data = data), 
     class = class)
}


# for internal use
# takes n-dimesional 'data' and gridlines to quantify each dimension
discretize_data = function(data, gridlines){
  
  # use gridlines for each dimension
  for(i in 1:ncol(data)){
    
    if(length(unique(gridlines))==0){
      discr = rep(1, nrow(data))
    }else{
      discr = rep(length(gridlines[[i]])+1, nrow(data))
      gridlines[[i]] = sort(gridlines[[i]])
      for(j in 1:length(gridlines[[i]])){ # determine discretization levels
        if(j == 1) {
          discr[data[,i] < gridlines[[i]][j]] = 1
        } else {
          discr[data[,i] < gridlines[[i]][j] & data[,i] >= gridlines[[i]][j-1]] = j
        }
      }
    }
    data[,i] = discr
  }
  
  return(data)
}
