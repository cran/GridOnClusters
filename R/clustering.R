# cluster.R
# 
# Jiandong Wang, Sajal Kumar, and Joe Song
#
# Created: Dec 12, 2025
#   Clustering related functions are extracted from 
#     discretize_jointly.R  

#' @title Cluster Multivariate Data
#' 
#' @description The function obtains clusters from data using the given
#'   number of clusters, which may be a range.
#'   
#' @param data input continuous multivariate data
#' @param k the number(s) of clusters 
#' @param method the method for clustering
#' @param noise adding jitter noise to the data or not
#'
#' @keywords internal
#' @export
cluster = function(data, k, method, noise)
{
  if(ncol(data) == 1) {
    modelName <- c("V")
  } else {
    # Gaussian balls of various radii
    modelName <- "VII"  
  }
  
  if(method == "Ball+BIC"){
    if(all(is.finite(k))){
      
      if(noise == TRUE){
        # use noisy Mclust
        mclust.res = Call_Noisy_Mclust(data, k)
        
      } else {
        
        mclust.res = mclust::Mclust(
          data, G = k, modelNames = modelName, 
          verbose = FALSE)
        
        # add additional information to Mclust (non noisy data)
        if(!is.null(mclust.res)){
          mclust.res$isNoisy = FALSE  
        }
        else{
          mclust.res = Call_Noisy_Mclust(data, k)
        }
      }
      
      cluster_info = list(
        clusters = as.vector(as.numeric(as.factor(mclust.res$classification))-1),
        data = as.matrix(data), method = "Ball+BIC")
      
    } else {
      
      range_k = c(2:8)
      
      # use noisy Mclust
      mclust.res = Call_Noisy_Mclust(data, range_k)
      
      # find the number of best clusters
      BIC_collect = mclust.res$BIC[,1]
      bestk = as.numeric(names(which.max(BIC_collect)))
      ogk = range_k
      while(bestk == max(range_k)){
        
        range_k = c((max(range_k)+1):min((max(range_k)*2), nrow(data)))
        ogk = c(ogk, range_k)
        
        mclust.res = Call_Noisy_Mclust(data, range_k)
        
        BIC_collect = c(BIC_collect, mclust.res$BIC[,1])  
        bestk = as.numeric(names(which.max(BIC_collect)))
        
      }
      
      # # let the user know that the user has reached the best k
      # message(paste0("Reached best k=",bestk,
      #                ". Re-clustering data using the best k."))
      
      # re-cluster using best k
      mclust.res = Call_Noisy_Mclust(data, bestk)
      
      cluster_info = list(
        clusters = as.vector(as.numeric(as.factor(mclust.res$classification))-1),
        data = as.matrix(data), method = method)
      
    }
  } else if(FALSE && method == "Ball+Tele+BIC"){
    # require(Telescope)
    # require(mclust)
    if(all(is.finite(k))){
      # tel = Telescope::telescope(data, TRUE)
      tele_tree = hc(data, modelName = modelName)
      tele_tree[,] = t(tele_tree) # t(my_tree)
      if(noise == TRUE){
        # use noisy Mclust
        mclust.res = Call_Noisy_Mclust(data, k) # , initialization = tele_tree)
      }else{
        mclust.res = mclust::Mclust(data, G = k, modelNames = modelName, 
                                    initialization = tele_tree, verbose = FALSE)
        
        # add additional information to Mclust (non noisy data)
        if(!is.null(mclust.res)){
          mclust.res$isNoisy = FALSE  
        }
      }
      
      cluster_info = list(
        clusters = as.vector(as.numeric(as.factor(mclust.res$classification))-1),
        data = as.matrix(data), method = "Ball+Tele+BIC")
      
    } else {
      
      range_k = c(2:8)
      
      # use noisy Mclust
      mclust.res = Call_Noisy_Mclust(data, range_k)
      
      # find the number of best clusters
      BIC_collect = mclust.res$BIC[,1]
      bestk = as.numeric(names(which.max(BIC_collect)))
      ogk = range_k
      while(bestk == max(range_k)){
        
        range_k = c((max(range_k)+1):min((max(range_k)*2), nrow(data)))
        ogk = c(ogk, range_k)
        
        mclust.res = Call_Noisy_Mclust(data, range_k)
        
        BIC_collect = c(BIC_collect, mclust.res$BIC[,1])  
        bestk = as.numeric(names(which.max(BIC_collect)))
        
      }
      
      # # let the user know that the user has reached the best k
      # message(paste0("Reached best k=",bestk,
      #                ". Re-clustering data using the best k."))
      
      # re-cluster using best k
      mclust.res = Call_Noisy_Mclust(data, bestk)
      
      cluster_info = list(
        clusters = as.vector(as.numeric(as.factor(mclust.res$classification))-1),
        data = as.matrix(data), method = method)
      
    }
    
  } 
  else if(method == "kmeans+silhouette"){
    # is k a single number or a range
    if(length(k) == 1 && is.finite(k)){
      
      # randomly generate centers
      centers = dqsample(which(!duplicated(data)), k)
      
      # get cluster information for 'k'
      kmeans.res = kmeans(data, centers = as.matrix(data[centers,]), 
                          iter.max = 100)
      
      # only keep cluster centers and labels
      cluster_info = list(#centers = as.matrix(kmeans.res$centers),
        clusters = as.vector(kmeans.res$cluster)-1,
        data = as.matrix(data), method = method)
      
    } else if(length(k) == 1 && !is.finite(k)) { # else if k is set to Inf
      
      data_dist = dist(data) # distance matrix for data
      
      # compute cluster info and silhouette score for the cluster range k, 
      # adaptively increasing k until optimal or each sample is placed in their 
      # own cluster
      
      # find unique number of samples
      unq = sum(!duplicated(data))
      if(unq > 8){
        range_k = c(2:8)  
      } else {
        range_k = c(2:unq)
      }
      
      # compute cluster info for original k range
      cluster_info = lapply(range_k, function(i){
        centers = dqsample(which(!duplicated(data)), i)
        data_clust = kmeans(data, centers = as.matrix(data[centers,]),
                            iter.max = 100)
        return(list(data_clust$cluster, data_clust$centers, 
                    mean(silhouette(data_clust$cluster, dist=data_dist)[,3])))
      })
      
      # find the max silhouette
      silhouette_scr = unlist(lapply(cluster_info, function(x){
        return(x[[3]])
      }), use.names = FALSE)
      max_silhouette = max(which.max(silhouette_scr)) # take the bigger k out of those with equal silhouette scores
      
      # while best k is equal to max k, keep computing scores
      bestk = range_k[max_silhouette]
      ogk = range_k
      while(bestk == max(range_k) && unq > 8){
        
        range_k = c((max(range_k)+1):min((max(range_k)*2), nrow(data)))
        ogk = c(ogk, range_k)
        
        # compute cluster info for original k range
        n_cluster_info = lapply(range_k, function(i){
          centers = dqsample(which(!duplicated(data)), i)
          data_clust = kmeans(data, centers = as.matrix(data[centers,]),
                              iter.max = 100)
          return(list(data_clust$cluster, data_clust$centers, 
                      mean(silhouette(data_clust$cluster, dist=data_dist)[,3])))
        })
        
        # append new cluster info to original
        cluster_info = c(cluster_info, n_cluster_info)
        
        # find the max silhouette
        silhouette_scr = unlist(lapply(cluster_info, function(x){
          return(x[[3]])
        }), use.names = FALSE)
        max_silhouette = max(which.max(silhouette_scr)) # take the bigger k out of those with equal silhouette scores
        
        # find best k
        bestk = ogk[max_silhouette]
      }
      
      # only keep cluster centers and labels for the 'k' with max silhouette
      cluster_info = list(#centers = as.matrix(cluster_info[[max_silhouette]][[2]]),
        clusters = as.vector(cluster_info[[max_silhouette]][[1]]-1),
        data = as.matrix(data), cluster_method = method)
      
    } else { # else k is a vector
      
      data_dist = dist(data) # distance matrix for data
      
      # compute cluster info and silhouette score for the cluster range k
      cluster_info = lapply(k, function(i){
        centers = dqsample(which(!duplicated(data)), i)
        data_clust = kmeans(data, centers = as.matrix(data[centers,]),
                            iter.max = 100)
        return(list(data_clust$cluster, data_clust$centers, 
                    mean(silhouette(data_clust$cluster, dist=data_dist)[,3])))
      })
      
      # find the max silhouette
      silhouette_scr = unlist(lapply(cluster_info, function(x){
        return(x[[3]])
      }), use.names = FALSE)
      max_silhouette = max(which.max(silhouette_scr)) # take the bigger k out of those with equal silhouette scores
      
      # only keep cluster centers and labels for the 'k' with max silhouette
      cluster_info = list(#centers = as.matrix(cluster_info[[max_silhouette]][[2]]),
        clusters = as.vector(cluster_info[[max_silhouette]][[1]]-1),
        data = as.matrix(data), cluster_method = method)
      
    }
  }
  else if(method == "PAM"){
    if(length(k) == 1 && is.finite(k)){
      # pre-cluster the data using partition around medoids (PAM)
      cluster_label = pam(x=data, diss = FALSE, metric = "euclidean", k = k)$clustering
    }
    else{
      # error
    }
    
    cluster_info = list(
      clusters = as.vector(cluster_label), data = as.matrix(data), method = method)
  }
  else if(FALSE && method == "telescope"){
    # require(Telescope)
    res <- NULL # Telescope::telescope(data,0)
    cluster_info = list(
      clusters = res$c_label, data = as.matrix(data), method = method)
  }
  return(cluster_info)
}

# for internal use
# takes n-dimensional 'data' and clusters it into 'k' groups using Mclust
# adding noise whenever necessary
Call_Noisy_Mclust = function(data, k){
  
  if(ncol(data) == 1) {
    modelName <- c("V")
  } else {
    # Gaussian balls of various radii
    modelName <- "VII"  
  }
  
  mclust.res = mclust::Mclust(data, G = k, modelNames = modelName, 
                              verbose = FALSE)
  
  # add additional information to Mclust (non noisy data)
  if(!is.null(mclust.res)){
    mclust.res$isNoisy = FALSE  
  }
  
  # Mclust can return NULL result for sparse discrete data
  amt = 0.5
  noisy_data = data # don't change the original input
  while(is.null(mclust.res)){
    
    # add a little bit of noise to enable clustering
    for(j in 1:ncol(noisy_data)){
      noisy_data[,j] = jitter(as.numeric(noisy_data[,j]), amount = amt)
    }
    
    mclust.res = mclust::Mclust(noisy_data, G = k, modelNames = modelName, 
                                verbose = FALSE)
    
    # add additional information to Mclust (noisy data)
    if(!is.null(mclust.res)){
      mclust.res$isNoisy = TRUE
    }
    
    # double the amount of noise for the next iteration (if any)
    amt = amt*2
    
  }
  
  if(mclust.res$isNoisy){
    warning(paste0("Mclust could not cluster 'data', jitter noise 
            (factor=1, amount=",(amt/2),") was added."))
  }
  
  return(mclust.res)
  
}

