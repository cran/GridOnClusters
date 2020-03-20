# discretize.jointly.R
#
# created by Sajal Kumar
# Copyright (c) NMSU Song lab

#' Discretize continuous multivariate data by a cluster-preserving grid
#'
#' Discretize continuous multivariate data using a grid that captures the joint distribution via
#' preserving clusters in the original data
#'
#' @importFrom GDAtools medoids
#' @importFrom cluster silhouette
#' @importFrom stats kmeans
#' @importFrom stats dist
#' @import Rcpp
#' @useDynLib GridOnClusters
#'
#' @param data a matrix containing two or more continuous variables.
#' Columns are variables, rows are observations.
#'
#' @param k either the number or range of clusters to be found on \code{data}.
#' The default is 2 to 10 clusters. If a range is specified, an optimal k in
#' the range is chosen to maximize the average silhouette width.
#' If \code{cluster_label} is specified, \code{k} is ignored.
#'
#' @param cluster_label a vector of user-specified cluster labels for each observation
#' in \code{data}. The user is free to choose any clustering.
#' If unspecified, k-means clustering is used by default.
#'
#' @return
#'
#' A list that contains two items:
#' \item{\code{D}}{a matrix that contains the discretized version of the original \code{data}.
#' The discretize values are one(1)-based.}
#'
#' \item{\code{grid}}{a list of vectors containing decision boundaries for each variable/dimension}
#'
#' @examples
#' # using a specified \code{k}
#' x = rnorm(100)
#' y = sin(x)
#' z = cos(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=5)$D
#'
#' # using a range of 'k'
#' x = rnorm(1000)
#' y = log1p(abs(x))
#' z = tan(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=c(3:10))$D
#'
#' # using an alternate different clustering scheme
#' library(cluster)
#' x = rnorm(1000)
#' y = log1p(abs(x))
#' z = sin(x)
#' data = cbind(x, y, z)
#'
#' # pre-cluster the data using partition around medoids (PAM)
#' cluster_label = pam(x=data, diss = FALSE, metric = "euclidean", k = 5)$clustering
#' discretized_data = discretize.jointly(data, cluster_label = cluster_label)$D
#'
#' @seealso
#'
#' See [`Ckmeans.1d.dp::Ckmeans.1d.dp()`] for discretizing a single continuous variable.
#'
#' @export
discretize.jointly = function(data, k=c(2:10), cluster_label=NULL){

  # check if data provided is a matrix
  if( !("matrix" %in% class(data)) && !("data.frame" %in% class(data))){
    stop("'data' must be a matrix or data.frame.")
  }

  # check if all columns in 'data' are numeric or integer
  dim_class = apply(data, 2, class)
  if(!all(dim_class == "numeric" | dim_class == "integer")){
    stop("All columns in 'data' should be numeric or integer.")
  }

  # 'data' should have atleast 10 points
  dim_data = dim(data)
  if(dim_data[1] < 10){
    stop("'data' should have atleast 10 observations,")
  }

  # 'cluster_label' should either be null or integers (or numeric) matching nrow(data)
  if(!is.null(cluster_label) && length(cluster_label) != nrow(data)){
    stop("'cluster_label' should either be null or a vector with nrow(data) elements.")
  }

  if(!is.null(cluster_label) && !class(cluster_label) %in% c("numeric","integer")){
    stop("'cluster_label' should be either null or a numeric/integer vector.")
  }

  # if no cluster labels are supplied, default to K-means
  if(is.null(cluster_label)){

    # is k a single number or a range
    if(length(k) == 1){

      # get cluster information for 'k'
      cluster_info = kmeans(data, centers = k)

      # only keep cluster centers and labels
      cluster_info = list(centers = as.matrix(cluster_info$centers),
                          clusters = as.vector(cluster_info$cluster)-1,
                          data = as.matrix(data))

    } else {

      data_dist = dist(data) # distance matrix for data

      # compute cluster info and silhouette score the cluster range k
      cluster_info = lapply(k, function(i){
        data_clust = kmeans(data, centers = i)
        return(list(data_clust$cluster, data_clust$centers, mean(silhouette(data_clust$cluster, dist=data_dist)[,3])))
      })

      # find the max silhouette
      silhouette_scr = unlist(lapply(cluster_info, function(x){
        return(x[[3]])
      }), use.names = FALSE)
      max_silhouette = max(which.max(silhouette_scr)) # take the bigger k out of those with equal silhouette scores

      # only keep cluster centers and labels for the 'k' with max silhouette
      cluster_info = list(centers = as.matrix(cluster_info[[max_silhouette]][[2]]),
                          clusters = as.vector(cluster_info[[max_silhouette]][[1]]-1),
                          data = as.matrix(data))
    }
  } else { # cluster labels are supplied

    cluster_label = as.numeric(as.factor(cluster_label)) # making labels consecutive

    # find medoids
    centers = data[medoids(D = dist(data), cl = cluster_label),,drop=FALSE]

    cluster_info = list(centers = as.matrix(centers),
                        clusters = cluster_label-1,
                        data = as.matrix(data))
  }

  # get grid lines
  grid_lines = findgrid(cluster_info, nrow(cluster_info$centers), nrow(data), ncol(data))

  # discretize data
  discr_data = discretize_data(data, grid_lines)

  return(list(D=discr_data, grid=grid_lines))
}

# for internal use
# takes n-dimesional 'data' and gridlines to quantify each dimension
discretize_data = function(data, gridlines){

  # use gridlines for each dimension
  for(i in 1:ncol(data)){

    discr = rep(length(gridlines[[i]])+1, nrow(data))
    gridlines[[i]] = sort(gridlines[[i]])
    for(j in 1:length(gridlines[[i]])){ # determine discretization levels
      if(j == 1) {
        discr[data[,i] < gridlines[[i]][j]] = 1
      } else {
        discr[data[,i] < gridlines[[i]][j] & data[,i] >= gridlines[[i]][j-1]] = j
      }
    }
    data[,i] = discr
  }

  return(data)
}
