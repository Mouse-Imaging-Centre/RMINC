# Voxel-wise ANOVA

Compute a sequential ANOVA at each voxel

## Usage

``` r
mincAnova(
  formula,
  data = NULL,
  subset = NULL,
  mask = NULL,
  maskval = NULL,
  parallel = NULL,
  cleanup = TRUE,
  conf_file = getOption("RMINC_BATCH_CONF")
)
```

## Arguments

- formula:

  The anova formula. The left-hand term consists of the MINC filenames
  over which to compute the models at every voxel.

- data:

  The dataframe which contains the model terms.

- subset:

  Subset definition.

- mask:

  Either a filename or a vector of values of the same length as the
  input files. ANOVA will only be computed inside the mask.

- maskval:

  the value in the mask used to select unmasked voxels, defaults to any
  positive intensity from 1-99999999 internally expanded to .5 -
  99999999.5. If a number is specified voxels with intensities within
  0.5 of the chosen value are considered selected.

- parallel:

  how many processors to run on (default=single processor). Specified as
  a two element vector, with the first element corresponding to the type
  of parallelization, and the second to the number of processors to use.
  For local running set the first element to "local" or "snowfall" for
  back-compatibility, anything else will be run with batchtools see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/pMincApply.md)
  Leaving this argument NULL runs sequentially.

- cleanup:

  Whether or not to remove parallelization files

- conf_file:

  A batchtools configuration file defaulting to
  `getOption("RMINC_BATCH_CONF")`

## Value

Returns an array with the F-statistic for each model specified by
formula with the following attributes:

- model design matrix

- filenames minc file names input

- dimensions dimensions of the statistics matrix

- dimnames names of the dimensions for the statistic matrix

- stat-type types of statistic used

- df degrees of freedom of each statistic

## Details

This function computes a sequential ANOVA over a set of files.

## See also

mincWriteVolume,mincFDR,mincMean, mincSd

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
# read the text file describing the dataset
gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
# run an ANOVA at each voxel
vs <- mincAnova(jacobians_fixed_2 ~ Sex, gf)
} # }
```
