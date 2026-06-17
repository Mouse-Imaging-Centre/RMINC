# False Discovery Rates

Takes the output of a minc modelling function and computes False
Discovery Rate thresholds.

## Usage

``` r
mincFDR(buffer, ...)

# S3 method for class 'mincSingleDim'
mincFDR(buffer, df, mask = NULL, method = "fdr", ...)

# S3 method for class 'mincLogLikRatio'
mincFDR(buffer, mask = NULL, ...)

# S3 method for class 'mincLmer'
mincFDR(buffer, mask = NULL, method = "fdr", ...)

# S3 method for class 'mincMultiDim'
mincFDR(
  buffer,
  columns = NULL,
  mask = NULL,
  df = NULL,
  method = "FDR",
  statType = NULL,
  ...
)
```

## Arguments

- buffer:

  The results of a mincLm type run.

- ...:

  extra parameters to pass to methods

- df:

  The degrees of freedom - normally this can be determined from the
  input object.

- mask:

  Either a filename or a numeric vector representing a mask only values
  inside the mask will be used to compute the threshold.

- method:

  The method used to compute the false discovery rate. Options are "FDR"
  and "pFDR".

- columns:

  A vector of column names. By default the threshold will be computed
  for all columns; with this argument the computation can be limited to
  a subset.

- statType:

  This should be either a "t","F","u","chisq" or "tlmer" depending upon
  the type of statistic being thresholded.

## Value

A object of type `mincQvals` with the same number of columns as the
input (or the subset specified by the columns argument to mincFDR). Each
column now contains the qvalues for each voxel. Areas outside the mask
(if a mask was specified) will be represented by a value of 1. The
result also has an attribute called "thresholds" which contains the 1,
5, 10, 15, and 20 percent false discovery rate thresholds.

## Details

This function uses the `qvalue` package to compute the False Discovery
Rate threshold for the results of a
[mincLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincLm.md)
computation. The False Discovery Rate represents the percentage of
results expected to be a false positive. Two implementations can be used
as specified by the method argument. "FDR" uses the implementation in
`p.adjust`, whereas "pFDR" is a version of the postivie False Discovery
Rate as found in John Storey's `qvalue` package. The main interface
functions are

- mincFDR.mincMultiDim The workhorse function, used to compute q-values
  and thresholds for sets of minc volumes

- mincFDR.logLikRatio Similar to above, but calculates thresholds by
  parametric bootstrap when possible

- mincFDR.mincSingleDim Used when
  [mincLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincLm.md)-like
  results are written out and read back in either to the same or another
  R session. In this case it loses it's `minMultiDim` class and must be
  converted back

- vertexFDR Used with results of a
  [vertexLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/vertexLm.md)-like
  command. Results are converted internally to resemble a mincMultiDim
  and processed as normal

## Methods (by class)

- `mincFDR(mincSingleDim)`: mincSingleDim

- `mincFDR(mincLogLikRatio)`: mincLogLikRatio

- `mincFDR(mincLmer)`: mincLmer

- `mincFDR(mincMultiDim)`: mincMultiDim

## See also

mincWriteVolume,mincLm,mincWilcoxon or mincTtest

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
# read the text file describing the dataset
gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
# compute the False Discovery Rate
qvals <- mincFDR(vs)
} # }
```
