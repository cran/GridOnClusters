# test-GridOnClusters.R
#
# tests discretize.jointly function under the perameter (cluster_method == "Ball+BIC" 
# and grid_method == "DP Compressed likelihood")
# Created by: Jiandong Wang, Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 16th August, 2022

library(testthat)
library(FunChisq)
library(GridOnClusters)
library(cluster)
library(dqrng)
library(mclust)

context("Testing Ball+BIC & DP approx likelihood 1-way")

test_that("Testing discretize.jointly (\"Ball+BIC\")(\"DP approx likelihood 1-way\")", {
   
   # test 1
   # y = f(x)
   # z = f(x)
   # k = constant
   
   cluster_method <- "Ball+BIC"
   grid_method <- "DP approx likelihood 1-way"
   
   dqset.seed(123)
   x = dqrnorm(100, mean=5, sd=1)
   y = sin(x)
   z = cos(x)
   data = cbind(x, y, z)
   discr = discretize.jointly(
      data, k=3, cluster_method=cluster_method,
      grid_method = grid_method, min_level = 1)
   
   # test marginal levels
   expect_equivalent(length(unique(discr$D[,1])), 4)
   expect_equivalent(length(unique(discr$D[,2])), 2)
   expect_equivalent(length(unique(discr$D[,3])), 3)
   
   # test marginal distribution
   expect_equivalent(table(discr$D[,1]), as.table(c(24, 39, 33, 4)))
   expect_equivalent(table(discr$D[,2]), as.table(c(46, 54)))
   expect_equivalent(table(discr$D[,3]), as.table(c(23, 43, 34)))
   
   # test 2d joint distributions
   dim12 = table(discr$D[,1], discr$D[,2])
   expect_equivalent(dim12, as.table(matrix(c(7, 17,
                                              39, 0,
                                              0, 33,
                                              0, 4),
                                            nrow=4, ncol=2, byrow = T)))
   
   dim13 = table(discr$D[,1], discr$D[,3])
   expect_equivalent(dim13, as.table(matrix(c(23, 1, 0,
                                              0, 39, 0,
                                              0, 0, 33,
                                              0, 3, 1),
                                            nrow=4, ncol=3, byrow = T)))
   
   dim23 = table(discr$D[,2], discr$D[,3])
   expect_equivalent(dim23, as.table(matrix(c(7, 39, 0,
                                              16, 4, 34),
                                            nrow=2, ncol=3, byrow = T)))
   
   # test ARI score
   #expect_equivalent(round(discr$csimilarity, digits = 3), 0.949)
   expect_equivalent(round(discr$csimilarity, digits = 3), 0.89)
   
   # test scale
   discr.scale = discretize.jointly(
      data, k=3, cluster_method=cluster_method,
      grid_method = grid_method, min_level = 1, scale = TRUE)
   
   # test marginal levels
   expect_equivalent(length(unique(discr.scale$D[,1])), 4)
   expect_equivalent(length(unique(discr.scale$D[,2])), 3)
   expect_equivalent(length(unique(discr.scale$D[,3])), 3)
   
   # test marginal distribution
   expect_equivalent(table(discr.scale$D[,1]), as.table(c(24, 39, 32, 5)))
   expect_equivalent(table(discr.scale$D[,2]), as.table(c(46, 47, 7)))
   expect_equivalent(table(discr.scale$D[,3]), as.table(c(23, 43, 34)))
   
   # test 2d joint distributions
   dim12 = table(discr.scale$D[,1], discr.scale$D[,2])
   expect_equivalent(dim12, as.table(matrix(c(7, 15, 2,
                                              39, 0, 0,
                                              0, 32, 0,
                                              0, 0, 5),
                                            nrow=4, ncol=3, byrow = T)))
   
   dim13 = table(discr.scale$D[,1], discr.scale$D[,3])
   expect_equivalent(dim13, as.table(matrix(c(23, 1, 0,
                                              0, 39, 0,
                                              0, 0, 32,
                                              0, 3, 2),
                                            nrow=4, ncol=3, byrow = T)))
   
   dim23 = table(discr.scale$D[,2], discr.scale$D[,3])
   expect_equivalent(dim23, as.table(matrix(c(7, 39, 0,
                                              15, 0, 32,
                                              1, 4, 2),
                                            nrow=3, ncol=3, byrow = T)))
   
   # test ARI score
   #expect_equivalent(round(discr.scale$csimilarity, digits = 3), 0.949)
   expect_equivalent(round(discr$csimilarity, digits = 3), 0.89)
   
   # test Scaled Upsilon score
   expect_equivalent(round(discr$upsilon_relt, digits = 3), 1.876)
   # test 2
   # y = f(x)
   # z = f(x)
   # k = variable (determined by silhouette)
   dqset.seed(321)
   x = dqrnorm(n = 100, mean=10, sd=2)
   y = log(x)
   z = tan(x)
   data = cbind(x, y, z)
   discr = discretize.jointly(data, k=c(3:10), cluster_method=cluster_method,
                              grid_method = grid_method, min_level = 1)
   
   # test marginal levels
   expect_equivalent(length(unique(discr$D[,1])), 6)
   expect_equivalent(length(unique(discr$D[,2])), 6)
   expect_equivalent(length(unique(discr$D[,3])), 3)
   
   # test marginal distribution
   expect_equivalent(table(discr$D[,1]), as.table(c(19, 19, 18, 12, 21, 11)))
   expect_equivalent(table(discr$D[,2]), as.table(c(19, 19, 18, 12, 21, 11)))
   expect_equivalent(table(discr$D[,3]), as.table(c(35, 28, 37)))
   
   # test 2d joint distributions
   dim12 = table(discr$D[,1], discr$D[,2])
   expect_equivalent(dim12, as.table(matrix(c(19,  0,  0,  0,  0,  0,
                                               0, 19,  0,  0,  0,  0,
                                               0,  0, 18,  0,  0,  0,
                                               0,  0,  0, 12,  0,  0,
                                               0,  0,  0,  0, 21,  0,
                                               0,  0,  0,  0,  0, 11),
                                            nrow=6, ncol=6, byrow = T)))
   dim13 = table(discr$D[,1], discr$D[,3])
   expect_equivalent(dim13, as.table(matrix(c(11,  2,  6,
                                               0, 19,  0,
                                               0,  0, 18,
                                               7,  0,  5,
                                              16,  5,  0,
                                               1,  2,  8),
                                            nrow=6, ncol=3, byrow = T)))
   dim23 = table(discr$D[,2], discr$D[,3])
   expect_equivalent(dim23, as.table(matrix(c(11,  2,  6,
                                              0, 19,  0,
                                              0,  0, 18,
                                              7,  0,  5,
                                              16,  5,  0,
                                              1,  2,  8),
                                            nrow=6, ncol=3, byrow = T)))
   
   # test ARI score
   #expect_equivalent(round(discr$csimilarity, digits = 3), 0.857)
   expect_equivalent(round(discr$csimilarity, digits = 3), 0.58)
   
   # test Scaled Upsilon score
   expect_equivalent(round(discr$upsilon_relt, digits = 3), 0.929)
   # test 3
   # y != f(x)
   # z = f(x, y)
   # k = variable (determined by silhouette)
   dqset.seed(1234)
   
   x = dqrexp(n=50, rate = 0.6)
   y = dqrnorm(50, mean=2, sd=0.5)
   z = sin(x) + cos(y)
   data = cbind(x, y, z)
   discr = discretize.jointly(
      data, k=c(3:10), min_level = 2, cluster_method=cluster_method,
      grid_method = grid_method)
   
   # test marginal levels
   expect_equivalent(length(unique(discr$D[,1])), 2)
   expect_equivalent(length(unique(discr$D[,2])), 2)
   expect_equivalent(length(unique(discr$D[,3])), 3)
   
   # test marginal distribution
   expect_equivalent(table(discr$D[,1]), as.table(c(27, 23)))
   expect_equivalent(table(discr$D[,2]), as.table(c(32, 18)))
   expect_equivalent(table(discr$D[,3]), as.table(c(15, 23, 12)))
   
   # test 2d joint distributions
   dim12 = table(discr$D[,1], discr$D[,2])
   expect_equivalent(
     dim12, as.table(matrix(c(19, 8,
                              13, 10),
                            nrow=2, ncol=2, byrow = T)))
   
   dim13 = table(discr$D[,1], discr$D[,3])
   expect_equivalent(
     dim13, as.table(matrix(c(5, 13, 9,
                              10, 10, 3),
                            nrow=2, ncol=3, byrow = T)))
   
   dim23 = table(discr$D[,2], discr$D[,3])
   expect_equivalent(
     dim23, as.table(matrix(c(2, 18, 12,
                              13, 5, 0),
                            nrow=2, ncol=3, byrow = T)))
   
   # test ARI score
   expect_equivalent(round(discr$csimilarity, digits = 3), 0.546)
   
   # test Scaled Upsilon score
   expect_equivalent(round(discr$upsilon_relt, digits = 3), 0.358)
   # test 4
   # y = f(x)
   # z = f(x)
   # k = fixed
   # using an alternate clustering strategy
   dqset.seed(2468)
   
   x = dqrnorm(n = 1000, mean = 10, sd = 2)
   y = sin(x)
   z = cos(y)
   data = cbind(x, y, z)
   # use PAM to cluster
   alt.cluster = pam(x = data, k = 5, diss = FALSE, metric = "euclidean", cluster.only = TRUE)
   
   discr = discretize.jointly(data = data, cluster_label = alt.cluster,
                              grid_method = grid_method, min_level = 1)
   
   # test marginal levels
   expect_equivalent(length(unique(discr$D[,1])), 5)
   expect_equivalent(length(unique(discr$D[,2])), 5)
   expect_equivalent(length(unique(discr$D[,3])), 4)
   
   # test marginal distribution
   expect_equivalent(table(discr$D[,1]), as.table(c(93, 195, 203, 320, 189)))
   expect_equivalent(table(discr$D[,2]), as.table(c(328, 25, 354, 149, 144)))
   expect_equivalent(table(discr$D[,3]), as.table(c(336, 222, 120, 322)))
   
   # test 2d joint distributions
   dim12 = table(discr$D[,1], discr$D[,2])
   expect_equivalent(dim12, as.table(matrix(c(12, 3, 41, 37, 0,
                                              0,  0,  0, 74, 121,
                                              0, 17,186, 0,  0, 
                                              315, 5, 0, 0,  0,
                                              1, 0, 127, 38, 23),
                                            nrow=5, ncol=5, byrow = T)))
   
   dim13 = table(discr$D[,1], discr$D[,3])
   expect_equivalent(dim13, as.table(matrix(c(6,   28,  20,  39,
                                              120, 46,  29,   0,
                                              0,    0,  34, 169,
                                              187,128,   5,   0,
                                              23,  20,  32, 114),
                                            nrow=5, ncol=4, byrow = T)))
   
   dim23 = table(discr$D[,2], discr$D[,3])
   expect_equivalent(dim23, as.table(matrix(c(194, 134,   0,   0,
                                              0,     0,  25,   0,
                                              0,     0,  34, 320,
                                              0,    86,  61,   2,
                                              142,   2,   0,   0),
                                            nrow=5, ncol=4, byrow = T)))
   
   # test ARI score
   expect_equivalent(round(discr$csimilarity, digits = 3), 0.61)
   
   # test Scaled Upsilon score
   expect_equivalent(round(discr$upsilon_relt, digits = 3), 0.753)
})
