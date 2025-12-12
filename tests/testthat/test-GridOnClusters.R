# test-GridOnClusters.R
#
# tests discretize.jointly function
# Created by: Jiandong Wang, Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 9th March, 2020

library(testthat)
library(FunChisq)
library(GridOnClusters)
library(cluster)
library(dqrng)

context("Testing kmeans+silhouette & Sort+split")

test_that("Testing discretize.jointly (\"kmeans+silhouette\")(\"Sort+split\")", {

  # test 1
  # y = f(x)
  # z = f(x)
  # k = constant

  cluster_method <- "kmeans+silhouette"
  grid_method <- "Sort+split"
  
  dqset.seed(123)
  x = dqrnorm(100, mean=5, sd=1)
  y = sin(x)
  z = cos(x)
  data = cbind(x, y, z)
  discr = discretize.jointly(
     data, k=3, cluster_method=cluster_method,
     grid_method = grid_method, min_level = 1)

  # test marginal levels
  expect_equivalent(length(unique(discr$D[,1])), 3)
  expect_equivalent(length(unique(discr$D[,2])), 2)
  expect_equivalent(length(unique(discr$D[,3])), 3)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(23, 45, 32)))
  expect_equivalent(table(discr$D[,2]), as.table(c(40, 60)))
  expect_equivalent(table(discr$D[,3]), as.table(c(26, 42, 32)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(5, 18,
                                             35,10,
                                             0, 32),
                                           nrow=3, ncol=2, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(22,  1,  0,
                                             4,  38,  3,
                                             0,   3, 29),
                                           nrow=3, ncol=3, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(9, 31, 0,
                                             17, 11, 32),
                                           nrow=2, ncol=3, byrow = T)))

  # test ARI score
  #expect_equivalent(round(discr$csimilarity, digits = 3), 1)
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.664)
  
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
  expect_equivalent(length(unique(discr$D[,1])), 9)
  expect_equivalent(length(unique(discr$D[,2])), 9)
  expect_equivalent(length(unique(discr$D[,3])), 6)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(9, 10, 39,  3,  1,  1,  5, 22, 10)))
  expect_equivalent(table(discr$D[,2]), as.table(c(9, 10, 39,  3,  1,  1,  5, 22, 10)))
  expect_equivalent(table(discr$D[,3]), as.table(c(1, 1, 7, 30, 59, 2)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(9, 0, 0, 0, 0, 0, 0, 0, 0,
                                             0,10, 0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 39, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0, 3, 0, 0, 0, 0, 0,
                                             0, 0, 0, 0, 1, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 1, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 5, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0,22, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0,10),
                                           nrow=9, ncol=9, byrow = T)))
  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(0,  0,  0,  1,  8,  0,
                                             0,  0,  2,  8,  0,  0,
                                             0,  0,  0,  2, 37,  0,
                                             0,  0,  0,  0,  1,  2,
                                             1,  0,  0,  0,  0,  0,
                                             0,  1,  0,  0,  0,  0,
                                             0,  0,  5,  0,  0,  0,
                                             0,  0,  0, 18,  4,  0,
                                             0,  0,  0,  1,  9,  0),
                                           nrow=9, ncol=6, byrow = T)))
  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(0,  0,  0,  1,  8,  0,
                                             0,  0,  2,  8,  0,  0,
                                             0,  0,  0,  2, 37,  0,
                                             0,  0,  0,  0,  1,  2,
                                             1,  0,  0,  0,  0,  0,
                                             0,  1,  0,  0,  0,  0,
                                             0,  0,  5,  0,  0,  0,
                                             0,  0,  0, 18,  4,  0,
                                             0,  0,  0,  1,  9,  0),
                                           nrow=9, ncol=6, byrow = T)))

  # test ARI score
  #expect_equivalent(round(discr$csimilarity, digits = 3), 0.968)
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.861)
  

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
  expect_equivalent(length(unique(discr$D[,1])), 3)
  expect_equivalent(length(unique(discr$D[,2])), 2)
  expect_equivalent(length(unique(discr$D[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(29, 17, 4)))
  expect_equivalent(table(discr$D[,2]), as.table(c(34, 16)))
  expect_equivalent(table(discr$D[,3]), as.table(c(9, 41)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(21, 8,
                                             12, 5,
                                             1, 3),
                                           nrow=3, ncol=2, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(1, 28,
                                             4, 13,
                                             4, 0),
                                           nrow=3, ncol=2, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(4, 30,
                                             5, 11),
                                           nrow=2, ncol=2, byrow = T)))

  # test ARI score
  #expect_equivalent(round(discr$csimilarity, digits = 3), 0.821)
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.534)
  
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
  expect_equivalent(length(unique(discr$D[,2])), 3)
  expect_equivalent(length(unique(discr$D[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(93, 196, 203, 319, 189)))
  expect_equivalent(table(discr$D[,2]), as.table(c(334, 451, 215)))
  expect_equivalent(table(discr$D[,3]), as.table(c(439, 561)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(13, 61, 19,
                                             0, 40, 156,
                                             5, 198, 0,
                                             315, 4, 0,
                                             1, 148, 40),
                                           nrow=5, ncol=3, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(20, 73,
                                             138, 58,
                                             0, 203,
                                             250, 69,
                                             31, 158),
                                           nrow=5, ncol=2, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(261, 73,
                                             0, 451,
                                             178, 37),
                                           nrow=3, ncol=2, byrow = T)))

  # test ARI score
  #expect_equivalent(round(discr$csimilarity, digits = 3), 0.881)
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.766)
  

})
