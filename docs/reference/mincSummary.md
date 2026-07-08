# Minc Voxel Summary Functions

Compute the mean, standard deviation, sum, or variance at every voxel
across a a set of MINC volumes. An optional grouping variable will split
the computation by group rather than performing it across all volumes as
is the default.

## Usage

``` r
mincSummary(
  filenames,
  grouping = NULL,
  mask = NULL,
  method = "mean",
  maskval = NULL
)

mincMean(filenames, grouping = NULL, mask = NULL, maskval = NULL)

mincVar(filenames, grouping = NULL, mask = NULL, maskval = NULL)

mincSum(filenames, grouping = NULL, mask = NULL, maskval = NULL)

mincSd(filenames, grouping = NULL, mask = NULL, maskval = NULL)

mincCorrelation(filenames, grouping, mask = NULL, maskval = NULL)
```

## Arguments

- filenames:

  Filenames of the MINC volumes across which to create the descriptive
  statistic.

- grouping:

  Optional grouping - contains same number of elements as filenames; the
  results will then have the descriptive statistic computed separately
  for each group, or in the case of method = "correlation" this is the
  variable to correlate against.

- mask:

  A mask specifying which voxels are to be included in the summary.

- method:

  the type of summarys statistic to calculate for each voxel

- maskval:

  the value in the mask used to select unmasked voxels, defaults to any
  positive intensity from 1-99999999 internally expanded to .5 -
  99999999.5. If a number is specified voxels with intensities within
  0.5 of the chosen value are considered selected.

## Value

The output will be a single vector containing as many elements as there
are voxels in the input files. If a grouping factor was specified then
the output will be a matrix consisiting of as many rows as there were
voxels in the files, and as many columns as there were groups.

## Functions

- `mincMean()`: mean

- `mincVar()`: Variance

- `mincSum()`: Sum

- `mincSd()`: Standard Deviation

- `mincCorrelation()`: Correlation

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv")
mm <- mincMean(gf$jacobians_0.2)
ms <- mincSd(gf$jacobians_0.2)
mv <- mincVar(gf$jacobians_0.2,gf$Strain)
ms2 <- mincSum(gf$jacobians_0.2,gf$Strain)
} # }
```
