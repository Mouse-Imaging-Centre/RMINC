# Write a MINC volume to file

Volume Export

## Usage

``` r
mincWriteVolume(buffer, ...)

# S3 method for class 'mincSingleDim'
mincWriteVolume(buffer, output.filename, clobber = NULL, ...)

# S3 method for class 'mincMultiDim'
mincWriteVolume(
  buffer,
  output.filename,
  column = 1,
  like.filename = NULL,
  clobber = NULL,
  ...
)

# Default S3 method
mincWriteVolume(buffer, output.filename, like.filename, clobber = NULL, ...)
```

## Arguments

- buffer:

  The data to be written to file. Usually the result of
  [mincLm](https://mouse-imaging-centre.github.io/RMINC/reference/mincLm.md)
  or some such command

- ...:

  additional arguments to pass to methods

- output.filename:

  The filename to which to write the data to

- clobber:

  Overwrite existing output file when set to TRUE, will not overwrite
  when set to FALSE and will prompt when NULL

- column:

  Optional name of the column of a multidimensional MINC object to write
  out. By default the first column is used

- like.filename:

  An existing MINC filename which has the same dimensions as the data to
  be written out. Normally this information is stored inside MINC data
  objects

## Value

A list with the parameters of the minc volume written

## Details

Writes a MINC volume to file

This function takes numeric data, usually the results computed from one
of the other mincFunctions, and writes it to file so that it can be
viewed or manipulated with the standard MINC tools

## Methods (by class)

- `mincWriteVolume(mincSingleDim)`: mincSingleDim

- `mincWriteVolume(mincMultiDim)`: mincMultiDim

- `mincWriteVolume(default)`: default

## See also

mincWriteVolume,mincLm,mincFDR,mincMean,mincSd

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
# read the text file describing the dataset
gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
# run a linear model relating the data in all voxels to Sex
vs <- mincLm(gf$jacobians_fixed_2 ~ Sex, gf)
# write the results to file
mincWriteVolume(vs, "Fstat.mnc", "F-statistic")
} # }
```
