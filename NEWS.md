---
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: 72
---

# NEWS

## Version 0.3.2  (major public release)
```         
  2025-12-12
```

1. [MAJOR PUBLIC RELEASE] Created version 0.3.2 from 0.3.1.2.
1. [MINOR CHANGE] Revised DESCRIPTION.
1. [MINOR CHANGE] Updated documentation for the exported gen_simdata() function.
1. [MINOR CHANGE] Updated testthat tests in test-Ball-approx_likelihood.R
1. [MINOR CHANGE] Updated vignette Examples.Rmd
1. [IMPROVEMENT] Deleted the Makevars file, resulting in replacing 
fixed CXX14 to default CXX17 C++ compiler.
1. [MINOR CHANGE] Cleaned up warning messages in C++ code due to signed / unsigned mismatch
1. Commented out the #pragma STDC FENV_ACCESS ON in Joint_Grid.cpp.
1. [MINOR CHANGE] Moved 'static const vector<double> *spx_c;' from Cluster.h
   to Cluster.cpp.
1. [MINOR CHANGE] Updated documentation of function discretize.jointly()
1. [MINOR CHANGE] Commented out some options to the cluster_method argument in 
discretize.jointly(). They are methods not yet ready for release.
1. [BUG FIX] Now discretize.jointly() works with on univariate (1-D) data
   without quitting.

## Version 0.3.1.2
```         
  2025-12-08
```

1.  [Bug fix] When cutting in decreasing order, the cutting positions (lines)
    were also sorted in decreasing order, causing only the first line to be used
    when constructing the table. The lines are now reordered in increasing order
    to prevent this issue.
    
```         
  2025-12-02
```

1.  [NEW FEATURE] Implemented bidirectional cutting—attempt cuts in both increasing 
    and decreasing order, then select the best cut.
2.  [NEW FEATURE] Introduced a new option for the grid_method parameter in the R 
    interface, enabling users to apply a bidirectional evaluation strategy 
    that tests cuts in both increasing and decreasing order and selects the 
    best result.
    
## Version 0.3.1.1
```         
  2025-11-28
```

1.  [NEW FEATURE] Introduced a new return value dimensional_cll that stores the
    best BIC scores across all dimensional cuts.

## Version 0.3.1

```         
  2025-11-04
```

1.  [IMPROVEMENT] The function now uses std::thread::hardware_concurrency() to detect
    available CPU cores on all supported platforms (including Linux).
2.  [IMPROVEMENT] Added automatic validation to handle invalid or excessive thread
    numbers and ensure at least one thread is used even when system
    detection fails.
3.  [MINOR CHANGE] Corresponding documentation for thread parameter (`nthread`)
    updated.

```         
  2025-11-07
```

1.  [IMPROVEMENT] Improved numerical consistency: adjusted C++ backend to ensure
    identical double-precision results across all platforms (macOS ARM,
    macOS Intel, Linux, and Windows).

## Version 0.3.0.3

```         
  2025-10-23
```

1.  [BUG FIX] Minor issue fixed in c++. (boundary incorrect for dp-binomial
    version)

```         
  2025-8-21
```

1.  [BUG FIX] Minor issue fixed in c++. (boundary incorrect for dp-recursion
    version)

## Version 0.3.0.2

```         
  2025-2-06
```

1.  [BUG FIX] Added robustness to the C++ implementation, ensuring it correctly
    handles data sets where all data points belong to a single cluster
    (i.e., all labels are identical).”

```         
  2025-2-04
```

1.  [NEW FEATURE] Introduced a new evaluation metric, mid_dist, which computes the
    Euclidean distance between points and the median of their respective
    blocks.
2.  [IMPROVEMENT] Optimized the C code to enhance reliability, improve execution
    speed, and increase readability.

## Version 0.3.0.1

```         
  2024-7-20
```

1.  [BUG FIX] Fixed a bug in the C++ implementation related to purity calculation.
2.  [NEW FEATURE] Added additional evaluation scores to the return value.
3.  [NEW FEATURE] Included discretization labels in the return value.
4.  [IMPROVEMENT] Improved the plot function to display column names as x-axis or
    y-axis labels, if available.

## Version 0.3.0

```         
  2024-4-23
```

1.  [NEW FEATURE] added "telescope" as a clustering method

```         
  2024-4-15
```

1.  [IMPROVEMENT] Added dqRNGkind("xoroshiro128+") before dqset.seed(...) in testthat
    files to avoid crash by the new version of dqrng

```         
  2024-3-16
```

1.  [IMPROVEMENT] Added a checking if mclust result is NULL which could lead the code
    crash when the data is very small and dense as mclust cannot work
    on such data. If mclust result id NULL, we will add some jitter
    noise to the data to let the code working correctly, which makes
    the code more robust.

## Version 0.2.0

```         
  2023-10-03
```

1.  [NEW FEATURE] Implements multithreads in C++ to speedup the process. User are
    allowed to input the number of thread they want to use, if it exceed
    the max number of thread that machine has, it will automatically set
    to the max number of thread provided by function
    hardware_concurrency().
2.  [MINOR CHANGE] Updated the test file as many parameters has been added.
3.  [BUG FIX] Fixed some bugs.

## Version 0.1.7

```         
  2023-09-02
```

1.  [NEW FEATURE] Added a new parameter "scale", the default value is false. If the
    data variance on each dimension differ to much, you are strongly
    suggested set it to true to obtain a better clustering result.Please
    noted the scaled data will only be used in the process of
    clustering.

## Version 0.1.6.1

```         
  2022-09-13
```

1.  [IMPROVEMENT] Added "const" keyword to parameters of some functions.
2.  [IMPROVEMENT] Added some type_cast to improve the robust.
3.  [IMPROVEMENT] Some minor structure improvement and optimization.

## Version 0.1.6

```         
  2022-08-27
```

1.  [NEW FEATURE] Added a new parameter "entropy", the default value is false. If user
    set it to true, then the score for the dynamic programming matrix
    will be calculated based on entropy instead of likelihood.
2.  [NEW FEATURE] entropy score: N(D~x~) \* log(N(D~x~)/N(D)) likelihood score:
    N(D~x~) \* N(D)\* log(N(D~x~)/N(D))

## Version 0.1.5

```         
  2022-08-25
```

1.  [MINOR CHANGE] Made minor changes to eliminate C/C++ warning messages in BIC.h,
    BIC.cpp, Cutting_Cluster_dp.cpp, Cutting_Cluster_dp_compressed.cpp,
    Joint_Grid.h, Joint_Grid.cpp.

2022-08-20

1.  [BUG FIX] In Example 9, Ruby's dataset, when grid_method = "Sort+Split", it
    will cut at 0; however, no line should be put if the left side equal
    to the right side. Check C code.

2022-08-17

1.  [IMPROVEMENT] Renamed the grid_method "Dp compressed likelihood" to "DP exact
    likelihood", "Dp compressed likelihood nlog" to "Dp approx
    likelihood".

2.  [IMPROVEMENT] Removed ULONG_MAX in "Cut_by_density.cpp" to fix compiling issues
    under Linux. Similarly, DBL_MAX and INT_MIN were removed in
    Cutting_Cluster_dp_compressed.cpp.

3.  [BUG FIX] "#include <cassert>" was added to "Joint_Grid.cpp" to fix a
    compiling error under Linux.

4.  [BUG FIX] Removed the usage of "cref()" in "Cutting_Cluster_dp.cpp" and
    "Cutting_Cluster_dp_compressed.cpp" to avoid compiling error on
    Linux.

5.  [NEW FEATURE] Adding test files for "DP Compressed likelihood nlogn" and "DP
    Compressed likelihood "in testthat.

6.  [NEW FEATURE] Added the perimeter "noise" for user to determined whether adding
    jitter noise to original data or not. The default value is FALSE
    which no noise will be added. Only will be used when cluster_method
    = "BALL+BIC".

7.  [MAJOR NEW FEATURE] Added a new grid method "DP Compressed likelihood nlogn", which is
    a speed up version of "DP Compressed likelihood" by using Divide and
    Conquer when filling the matrix, but no optimal will be generated.
    The "DP Compressed likelihood nlogn" is set to default for the
    parameter "grid_method".

8.  [IMPROVEMENT] Optimized the structure of the cluster class in c code to speed up
    and make the code more robust, and rewrite the "Sort+Split" part
    correspondingly.

9.  [NEW FEATURE] Added cluster_method and grid_method to the result, if the cluster
    label is given by user, the cluster_method will be set to "user
    supply".

10. [IMPROVEMENT] Modify the plot function to print both methods in title.

## Version 0.1.4

2022-08-07

1.  [BUG FIX] Replaced ULONG_MAX by std::numeric_limits<unsigned long>::max() in
    "Cut_by_density.cpp" to fix compiling issues under Linux. Header
    file <limits> is included. Similarly, DBL_MAX and INT_MIN were
    replaced in Cutting_Cluster_dp_compressed.cpp.

2.  [BUG FIX] Included header file <functional> to fix compiling error on Linux
    regarding cref() in "Cutting_Cluster_dp.cpp" and
    "Cutting_Cluster_dp_compressed.cpp".

3.  [BUG FIX] "#include <cassert>" was added to "Joint_Grid.cpp" to fix a
    compiling error under Linux.

2022-07-27

4.  Version 0.1.4 created from version 0.1.3.

## Version 0.1.3

2022-06-30

1.  Version 0.1.3 Created from version 0.1.2.

1.  [NEW FEATURE] Used dynamic programming to maximize the likelihood of categorical
    distribution, to improve quality of discretization

1.  [NEW FEATURE] Used BIC on the categorical likelihood to select number of discrete
    levels for each variable.

1.  [IMPROVEMENT] Compression: avoid trying discretization within consecutive points
    of the same cluster label to greatly reduce runtime of dynamic
    programming.
    
1. [BUG FIX] The Examples.Rmd vignette does not compile on MacOS.

1. [BUG FIX] Five test cases failed.


## Version 0.1.2

2022-03-17

1.  [NEW FEATURE] Added an additional parameter cutoff, default value is 0.05. This
    parameter is only used for the "DP" and "DP Compressed". After
    calculating the line for each dimension, we will calculate the
    p-value. If the p-value is larger than the cutoff, no line will be
    kept for that dimension. Otherwise, store those lines.

## Version 0.1.1

2022-03-10

1.  [NEW FEATURE] Added new cluster methods "DP" and "DP Compressed". The new method
    using dynamic programming to find the optimal solution to split the
    label into many zone. The maximum number of zone is set to 2\*K by
    hand currently.

2.  [NEW FEATURE] Implemented a new approach of the "DP" method. Improved by
    compressed the label, the size of the dynamic programming table can
    be decreased dramatically.

3.  [IMPROVEMENT] Before returning the bins on each dimension, we will check the upsilon
    static to make sure the cutting is reasonable. Currently the cut off
    has been set to the median of the PDF. Any result with a upsilon
    static that smaller than the median of the PDF will be drop and no
    cut will be put on that dimension.

## Version 0.1.0

2022-01-25

1.  When cluster_method = "Ball+BIC", if Mclust discretization was
    unsuccessful, discretize.jointly will now add jitter noise (of
    increasing 'amount') to 'data'.

2022-01-18

1.  Updated CITATION, REFERENDCES.bib, and README.md.

2022-01-17

1.  Fixed a bug in computing weights for MultiChannel.WUC when cluster
    IDs are not consecutive numbers after Ball clustering by Mclust
    (discretize_jointly.R:198-204)
2.  Fixed function cluster(data, k, method) where the method was not
    included into the return object.
3.  Examples in the vignettes have been updated.

2022-01-11

1.  Added k=Inf option, that incrementally increases the number of
    clusters to choose until either an increment is not required or each
    data point is its own cluster.
2.  Added an example showing the use of k=Inf option.

2021-10-26

1.  Added an additional parameter 'cluster_method', could be "Ball+BIC",
    "kmeans+silhouette" or "PAM",
2.  Added clustering method Mclust
3.  Changed the default cluster method from Kmeans to Mclust, which is
    "Ball+BIC"
4.  Added an additional parameter 'grid_method', could be "Sort+split"
    or "MultiChannel.WUC", "MultiChannel.WUC" is a discretize method
    Multichannel.MUC from Ckmeans.1d.dp. The default option is still
    "Sort+split" which is the same as previous version. As the
    "MultiChannel.WUC" option is experimental, it is not recommended.
5.  Added general plotting function plot(), old plot function still
    available
6.  Changed the default value of min_level from 2 to 1. We expanded
    min_level argument to provide a vector of integers specifying the
    minimum level for each dimension.

## Version 0.0.8.1

2021-10-06

1.  Removed the requirement for the number of points in the data

## Version 0.0.8

2020-09-13

1.  Added function plotGOCpatterns to plot the continuous data along
    with the cluster preserving grid.
2.  Created a manual for the plotGOCpatterns() function.
3.  Updated the code for the Examples vignette to use plotGOCpatterns.

2020-08-10

1.  Created version 0.0.8 from 0.0.7.
2.  Added an additional parameter 'min_level' to denote the minimum
    number of discretization levels required for each dimension.
3.  Updated the manual of discrete.jointly() function.
4.  Added an entry in reference and citation.
5.  Updated README.md with badges.

## Version 0.0.7

2020-04-03

1.  Tidied up code for the Examples vignette.
2.  Updated the manual of discrete.jointly() function.
3.  Made minor editorial changes in DESCRIPTION and README.md.
4.  Resolved signed/unsigned mismatches.

2020-03-31

1.  Created version 0.0.7 from 0.0.6.
2.  Fixed memory leak in Clusters.cpp when calculating median.

## Version 0.0.6

2020-03-26

1.  Created version 0.0.6 from 0.0.5.
2.  Rewrote Prep_Index() to work in between two consecutive points,
    rather than on top of a single point.
3.  Using distance() in Prep_Index() to calculate the distance for two
    iterators.
4.  Using "ceil" in Binary_Index_Searching() to consider even/odd cases
    when determining grid lines.
5.  Fixed potential overflow issues.

## Version 0.0.5 (not released to the public)

2020-03-25

1.  Fixed a bug in the prep_index() function.
2.  Fixed prep_index() (lines 120 and 125) such that grid lines are put
    at the midpoint between two conseuctive points, instead of on one of
    the points.
3.  Updated vignette. Example 2 seems always correct now.

2020-03-24

1.  Created version 0.0.5 from 0.0.4.
2.  Function discretize.jointly() now returns cluster labels of each
    observation and a similarity score (ARI) between the joint
    discretization and the cluster labels of each observation.
3.  The class Cluster has a new constructor that takes cluster labels
    and the input data to compute median for each dimension.
4.  Find_grid() is now based on median.
5.  Using 'dqrng' in test cases to avoid RNG issue in testing.
6.  Rewrote multiple functions in Joint_Grid.cpp to avoid push_back().
7.  New visualization code in vignette now shows cluster labels for each
    observation.

## Version 0.0.4 (not released to the public)

2020-03-20

1.  Created version 0.0.4 from 0.0.3.
2.  Fixed typos in DESCRIPTION and README files.

## Version 0.0.3

2020-03-17

1.  Created version 0.0.3 from 0.0.2. Package renamed to GridOnClusters
2.  Function joint.grid.discretize.R() renamed to discretize.jointly()
3.  Return values of function discretize.jointly() changed to include
    both the discretized data and the grid
4.  Manual for discretize.jointly() updated.
5.  Line 104, 105 in Joint_Grid.cpp commented out
6.  Rewrote Find_Grid() to avoid push_back() in Joint_Grid.cpp
7.  Created a vignette to include examples.

## Version 0.0.2 (not released to the public)

2020-03-14

1.  Created the initial version 0.0.1. Package renamed to QNJGD

## Version 0.0.1 (not released to the public)

2020-03-09

1.  Created the initial version 0.0.1. Package named JointGridDiscr
