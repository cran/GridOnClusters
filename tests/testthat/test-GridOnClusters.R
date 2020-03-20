# test-GridOnClusters.R
#
# tests discretize.jointly function
# Created by: Jiandong Wang, Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 9th March, 2020

library(testthat)
library(FunChisq)
library(GridOnClusters)
library(cluster)

context("Testing GridOnClusters")

test_that("Testing discretize.jointly", {

  # test 1
  # y = f(x)
  # z = f(x)
  # k = constant
  set.seed(123)

  x = rnorm(100, mean=5, sd=1)
  y = sin(x)
  z = cos(x)

  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=3)$D

  # test marginal levels
  expect_equivalent(length(unique(discr[,1])), 3)
  expect_equivalent(length(unique(discr[,2])), 3)
  expect_equivalent(length(unique(discr[,3])), 3)

  # test marginal distribution
  expect_equivalent(table(discr[,1]), as.table(c(24, 49, 27)))
  expect_equivalent(table(discr[,2]), as.table(c(56, 18, 26)))
  expect_equivalent(table(discr[,3]), as.table(c(24, 48, 28)))

  # test functional dependency
  dim12 = table(discr[,1], discr[,2])
  expect_equivalent(dim12, as.table(matrix(c(12, 8, 4,
                                             44, 5, 0,
                                             0, 5, 22),
                                           nrow=3, ncol=3, byrow = T)))

  stat12 = fun.chisq.test(dim12)
  expect_equivalent(round(stat12$statistic, digits = 2), 80.54)
  expect_equivalent(signif(stat12$p.value, digits = 3), 1.34e-16)

  dim13 = table(discr[,1], discr[,3])
  expect_equivalent(dim13, as.table(matrix(c(24, 0, 0,
                                             0, 46, 3,
                                             0, 2, 25),
                                           nrow=3, ncol=3, byrow = T)))

  stat13 = fun.chisq.test(dim13)
  expect_equivalent(round(stat13$statistic, digits = 2), 162.07)
  expect_equivalent(signif(stat13$p.value, digits = 3), 5.26e-34)


  # test 2
  # y = f(x)
  # z = f(x)
  # k = variable (determined by silhouette)
  set.seed(321)

  x = rchisq(n = 50, df = 4)
  y = log(x)
  z = tan(x)

  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=c(3:10))$D

  # test marginal levels
  expect_equivalent(length(unique(discr[,1])), 7)
  expect_equivalent(length(unique(discr[,2])), 6)
  expect_equivalent(length(unique(discr[,3])), 5)

  # test marginal distribution
  expect_equivalent(table(discr[,1]), as.table(c(6, 2, 2, 12, 13, 9, 6)))
  expect_equivalent(table(discr[,2]), as.table(c(6, 2, 2, 26, 8, 6)))
  expect_equivalent(table(discr[,3]), as.table(c(2, 37, 7, 3, 1)))

  # test functional dependency
  dim12 = table(discr[,1], discr[,2])
  expect_equivalent(dim12, as.table(matrix(c(6, 0, 0, 0, 0, 0,
                                             0, 2, 0, 0, 0, 0,
                                             0, 0, 2, 0, 0, 0,
                                             0, 0, 0, 12, 0, 0,
                                             0, 0, 0, 13, 0, 0,
                                             0, 0, 0, 1, 8, 0,
                                             0, 0, 0, 0, 0, 6),
                                           nrow=7, ncol=6, byrow = T)))

  stat12 = fun.chisq.test(dim12)
  expect_equivalent(round(stat12$statistic, digits = 2), 190.93)
  expect_equivalent(signif(stat12$p.value, digits = 3), 2.43e-25)

  dim13 = table(discr[,1], discr[,3])
  expect_equivalent(dim13, as.table(matrix(c(0, 3, 3, 0, 0,
                                             0, 0, 0, 1, 1,
                                             2, 0, 0, 0, 0,
                                             0, 12, 0, 0, 0,
                                             0, 9, 4, 0, 0,
                                             0, 7, 0, 2, 0,
                                             0, 6, 0, 0, 0),
                                           nrow=7, ncol=5, byrow = T)))

  stat13 = fun.chisq.test(dim13)
  expect_equivalent(round(stat13$statistic, digits = 2), 43.55)
  expect_equivalent(signif(stat13$p.value, digits = 3), 0.0086)


  # test 3
  # y != f(x)
  # z = f(x, y)
  # k = variable (determined by silhouette)
  set.seed(1234)

  x = rchisq(n = 50, df = 4)
  y = rnorm(50, mean=2, sd=0.5)
  z = sin(x) + cos(y)

  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=c(3:10))$D

  # test marginal levels
  expect_equivalent(length(unique(discr[,1])), 4)
  expect_equivalent(length(unique(discr[,2])), 2)
  expect_equivalent(length(unique(discr[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr[,1]), as.table(c(19, 22, 8, 1)))
  expect_equivalent(table(discr[,2]), as.table(c(21, 29)))
  expect_equivalent(table(discr[,3]), as.table(c(18, 32)))

  # test functional dependency
  dim12 = table(discr[,1], discr[,2])
  expect_equivalent(dim12, as.table(matrix(c(8, 11,
                                             6, 16,
                                             6, 2,
                                             1, 0),
                                           nrow=4, ncol=2, byrow = T)))

  stat12 = fun.chisq.test(dim12)
  expect_equivalent(round(stat12$statistic, digits = 2), 6.74)
  expect_equivalent(signif(stat12$p.value, digits = 3), 0.0807)

  dim13 = table(discr[,1], discr[,3])
  expect_equivalent(dim13, as.table(matrix(c(1, 18,
                                             13, 9,
                                             4, 4,
                                             0, 1),
                                           nrow=4, ncol=2, byrow = T)))

  stat13 = fun.chisq.test(dim13)
  expect_equivalent(round(stat13$statistic, digits = 2), 13.02)
  expect_equivalent(signif(stat13$p.value, digits = 3), 0.0046)

  dim23 = table(discr[,2], discr[,3])
  expect_equivalent(dim23, as.table(matrix(c(4, 17,
                                             14, 15),
                                           nrow=2, ncol=2, byrow = T)))

  stat23 = fun.chisq.test(dim23)
  expect_equivalent(round(stat23$statistic, digits = 2), 4.16)
  expect_equivalent(signif(stat23$p.value, digits = 3), 0.0413)


  # test 4
  # y = f(x)
  # z = f(y)
  # k = variable (determined by silhouette)
  set.seed(4321)

  x = rchisq(n = 100, df = 5)
  y = sin(x)
  z = cos(y)

  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=c(3:10))$D

  # test marginal levels
  expect_equivalent(length(unique(discr[,1])), 7)
  expect_equivalent(length(unique(discr[,2])), 3)
  expect_equivalent(length(unique(discr[,3])), 3)

  # test marginal distribution
  expect_equivalent(table(discr[,1]), as.table(c(19, 14, 25, 13, 15, 12, 2)))
  expect_equivalent(table(discr[,2]), as.table(c(29, 39, 32)))
  expect_equivalent(table(discr[,3]), as.table(c(56, 25, 19)))

  # test functional dependency
  dim12 = table(discr[,1], discr[,2])
  expect_equivalent(dim12, as.table(matrix(c(0, 3, 16,
                                             0, 14, 0,
                                             23, 2, 0,
                                             0, 12, 1,
                                             0, 0, 15,
                                             6, 6, 0,
                                             0, 2, 0),
                                           nrow=7, ncol=3, byrow = T)))

  stat12 = fun.chisq.test(dim12)
  expect_equivalent(round(stat12$statistic, digits = 2), 148.68)
  expect_equivalent(signif(stat12$p.value, digits = 3), 1.05e-25)

  dim23 = table(discr[,2], discr[,3])
  expect_equivalent(dim23, as.table(matrix(c(28, 1, 0,
                                             0, 20, 19,
                                             28, 4, 0),
                                           nrow=3, ncol=3, byrow = T)))

  stat23 = fun.chisq.test(dim23)
  expect_equivalent(round(stat23$statistic, digits = 2), 91.09)
  expect_equivalent(signif(stat23$p.value, digits = 3), 7.74e-19)


  # test 5
  # y = f(x)
  # z = f(x)
  # k = fixed
  # using an alternate clustering strategy
  set.seed(2468)

  x = rnorm(n = 1000, mean = 10, sd = 2)
  y = sin(x)
  z = cos(y)

  data = cbind(x, y, z)

  # use PAM to cluster
  alt.cluster = pam(x = data, k = 5, diss = FALSE, metric = "euclidean", cluster.only = TRUE)

  discr = discretize.jointly(data = data, cluster_label = alt.cluster)$D

  # test marginal levels
  expect_equivalent(length(unique(discr[,1])), 5)
  expect_equivalent(length(unique(discr[,2])), 4)
  expect_equivalent(length(unique(discr[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr[,1]), as.table(c(155, 219, 219, 267, 140)))
  expect_equivalent(table(discr[,2]), as.table(c(275, 222, 314, 189)))
  expect_equivalent(table(discr[,3]), as.table(c(484, 516)))

  # test functional dependency
  dim12 = table(discr[,1], discr[,2])
  expect_equivalent(dim12, as.table(matrix(c(17, 8, 53, 77,
                                             0, 0, 136, 83,
                                             35, 141, 43, 0,
                                             223, 44, 0, 0,
                                             0, 29, 82, 29),
                                           nrow=5, ncol=4, byrow = T)))

  stat12 = fun.chisq.test(dim12)
  expect_equivalent(round(stat12$statistic, digits = 2), 1094.8)
  expect_equivalent(signif(stat12$p.value, digits = 3), 7.63e-227)

  dim13 = table(discr[,1], discr[,3])
  expect_equivalent(dim13, as.table(matrix(c(97, 58,
                                             109, 110,
                                             29, 190,
                                             211, 56,
                                             38, 102),
                                           nrow=5, ncol=2, byrow = T)))

  stat13 = fun.chisq.test(dim13)
  expect_equivalent(round(stat13$statistic, digits = 2), 246.39)
  expect_equivalent(signif(stat13$p.value, digits = 3), 3.9e-52)

})
