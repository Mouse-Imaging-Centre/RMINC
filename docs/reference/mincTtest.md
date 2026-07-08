# Minc T-test

Perform an unpaired,unequal variance t-test across a set of minc volumes

## Usage

``` r
mincTtest(filenames, grouping, mask = NULL, maskval = NULL)
```

## Arguments

- filenames:

  Filenames of the MINC volumes across which to run the t-test

- grouping:

  Contains same number of elements as filenames; must contain exactly
  two groups with which to compare means

- mask:

  A mask specifying which voxels are to be included in the test

- maskval:

  The value with which to mask the data (data will masked +/- 0.5 around
  this value

## Value

The output will be a single vector containing as many elements as there
are voxels in the input files, with that voxel's t-statistic

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv")
mtt <- mincTtest(gf$jacobians_0.2,gf$Strain)
} # }
```
