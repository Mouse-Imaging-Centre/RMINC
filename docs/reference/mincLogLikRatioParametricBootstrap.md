# computes parametric bootstrap on mincLogLikRatio output

The Log Likelihood Ratio tests closely approximates a Chi-squared
distribution when the number of groups (i.e. individual subjects in a
longitudinal study) is large (\>50), but can be anticonservative when
small. A parametric bootstrap test, in which data is randomly simulated
from the null model and then fit with both models, can give the correct
p-value. Here we compute the parametric boostrap on a small number of
randomly chosen voxels to get a sense of biased the estimated p-values
from the log likelihood ratio test really were.

## Usage

``` r
mincLogLikRatioParametricBootstrap(
  logLikOutput,
  selection = "random",
  nsims = 500,
  nvoxels = 50
)
```

## Arguments

- logLikOutput:

  the output from mincLogLikRatio

- selection:

  the algorithm for randomly chosing voxels. Only "random" works for
  now.

- nsims:

  the number of simulations to run per voxel

- nvoxels:

  the number of voxels to run the parametric bootstrap on

## Value

a matrix containing the chi-square p-values and the bootstrapped
p-values
