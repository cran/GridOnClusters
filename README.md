The 'GridOnClusters' R package
===============================

### Overview

The package offers a method to discretize multivariate continuous data using a grid that captures the joint distribution via preserving clusters in the original data. Joint grid discretization is applicable as a data transformation step before using other methods to infer association, function, or causality without assuming a parametric model.

### When to use the package

Most available discretization methods process one variable at a time, such as ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp). If discretizing each variable independently misses patterns arising from the joint distribution of multiple involved variables, one may benefit from using the joint discretization method in this package.

### To download and install the package

```{r}
install.packages("GridOnClusters")
```

### Examples

See the Examples vignette of the package.
