# GridOnClusters/tools.R
#
# Created  by Jiandong Wang
# Copyright (c) NMSU Song lab

#' @importFrom stats pchisq rnorm
#' @importFrom cluster pam
#' @importFrom mclust hc
#' 
#' @title Generate Simulated Data
#' @param cord data matrix that records the index for
#'   each cluster on each dimension
#' @param sim_table a matrix
#' @param noise a numeric value to specify noise level
#' @param plot a logical to turn on or off plotting
#' @keywords internal
#' @export
gen_simdata = function(cord, sim_table, noise=0.3, plot = FALSE){
  
   # n is the number of points of the sample
   n <- sum(sim_table)
   
   rsums = rowSums(sim_table)
   csums = colSums(sim_table)
   ndim = ncol(cord)
   # ndim does not include the Y axis, so total dimension will be ndim+1
   dist = 5- 5*noise
   data = matrix(nrow = n, ncol = ndim+1)
   
   label = c()
   cur_l = 1
   y = c()
   x = c()
   
   start_indices <- cumsum(c(1, rsums[-length(rsums)])) 
   end_indices <- cumsum(rsums)
   if(end_indices[nrow(sim_table)] != n){
      print("Error: end_indices does not match data size(n)!")
      return(-1)
   }
   
   if(ndim>=2){
      x_dim = c()
      for (i in 1:nrow(sim_table)) {
         x_temp = c()
         y_temp = c()
         
         for (j in 1:ndim) {
            dim.index = cord[i,j]
            if(start_indices[i] <= end_indices[i]){
               data[start_indices[i]:end_indices[i], j] = rnorm(rsums[i], dim.index*dist)
            }
         }
         
         for(j in 1:ncol(sim_table)){
            y_temp = c(y_temp, rnorm(sim_table[i,j], j*dist))
            label = c(label, rep(cur_l, sim_table[i,j]))
            cur_l = cur_l+1
         }
         y = c(y, y_temp)
      }
      data[,ndim+1] = y
   }
   else{
      for (i in 1:nrow(sim_table)) {
         x_sub = rnorm(rsums[i], dist*i)
         x = c(x, x_sub)
         for (j in 1:ncol(sim_table)) {
            label = c(label, rep(cur_l,sim_table[i,j]))
            y_sub = rnorm(sim_table[i,j],dist*j)
            y = c(y, y_sub)
            cur_l = cur_l+1
         }
      }
      # why max(y)+3-y? what this works?
      # y = max(y)+3-y
      if(plot){
         plot(x,y)
      }
      data = cbind(x,y)
   }
   
   res = structure(list(sim_data=data, label = as.numeric(as.factor(label))))
   return(res)
}

# NO @export
hashed_upsilon = function(points, log.p=FALSE){
   if((is.matrix(points) || is.data.frame(points)) &&
      nrow(points) > 1 && ncol(points) > 1) {
      res <- upsilon_c(as.matrix(points), log.p)
   } else {
      res <- list(0, 0, 1, 0)
   }
   
   max_upsilon <- (prod(apply(points, 2, function(col) length(unique(col))))*(nrow(points)))/(2^(ncol(points)))
   urelt <- res[[1]] / max_upsilon
   prelt <- pchisq(urelt, df = 1, lower.tail = FALSE, log.p = TRUE)
   res$stat_relt <- urelt
   res$p.value_relt <- prelt
   res$tsp <- prod(apply(points, 2, function(col) length(unique(col))))
   res$tsm <- mean(apply(points, 2, function(col) length(unique(col))))
   
   full_max <- 2 * sqrt(max_upsilon) -1
   urelt <- (res[[1]] / nrow(points)) / full_max
   
   names(res) <- c("statistic", "parameter", "p.value", "estimate", "stat_relt", "p.value_relt", "tsp", "tsm")
   names(res$statistic) <- "Upsilon"
   names(res$parameter) <- "df"
   names(res$estimate) <- "Effect size"
   names(res$p.value) <- ifelse(log.p, "log P", "P")
   names(res$stat_relt) <- "Upsilon relative"
   names(res$p.value_relt) <- ifelse(log.p, "log P relative", "P relative")
   names(res$tsp) <- "table_size_production"
   names(res$tsm) <- "table_size_mean"

   level <- res[[1]] / max_upsilon
   res$upsilon_level <- level
   return(res)
}

# NO @export
median_distance = function(data, label){
   dist.temp <- list()
   dist.temp$distance_to_median <- median_distance_c(data, label)
   dist.temp$mean_dist_to_median <- mean(unlist(dist.temp$distance))
   return(dist.temp)
}

# NO @export
mean_distance = function(data, label){
   dist.temp <- list()
   dist.temp$distance_to_mean <- mean_distance_c(data, label)
   dist.temp$mean_dist_to_mean <- mean(unlist(dist.temp$distance))
   return(dist.temp)
}

# NO @export
fill_columns <- function(mat) {
   existing_cols <- as.numeric(colnames(mat))
   full_range <- seq(min(existing_cols), max(existing_cols))  # Get full range of column indices
   missing_cols <- setdiff(full_range, existing_cols)  # Find missing columns
    
   if (length(missing_cols) > 0) {
      zero_mat <- matrix(0, nrow = nrow(mat), ncol = length(missing_cols))
      colnames(zero_mat) <- as.character(missing_cols)
      
      mat <- cbind(mat, zero_mat)
      mat <- mat[, order(as.numeric(colnames(mat)))]
   }
   return(mat)
}
