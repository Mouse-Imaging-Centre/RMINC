# Minc Wilcoxon

Perform a Mann-Whitney U test between a set of minc volumes.

## Usage

``` r
mincWilcoxon(filenames, grouping, mask = NULL, maskval = NULL)
```

## Arguments

- filenames:

  Filenames of the MINC volumes across which to run the test

- grouping:

  Contains same number of elements as filenames; must contain exactly
  two groups.

- mask:

  A mask specifying which voxels are to be included in the test

- maskval:

  The value with which to mask the data (data will masked +/- 0.5 around
  this value

## Value

The output will be a single vector containing as many elements as there
are voxels in the input files, with that voxel's U value (lower one).
The number of observations in each sample is also saved as an
attribute(m and n) so the result can be passed into mincFDR.

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv")
mw <- mincWilcoxon(gf$jacobians_0.2,gf $Strain)
} # }
```
