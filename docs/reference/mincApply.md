# Apply arbitrary R function at every voxel

Can execute any R function at every voxel for a set of MINC files

## Usage

``` r
mincApply(
  filenames,
  function.string,
  mask = NULL,
  maskval = NULL,
  reduce = FALSE
)
```

## Arguments

- filenames:

  The MINC files over which to apply the function. Have to be the same
  sampling.

- function.string:

  The function which to apply. Can only take a single argument, which
  has to be 'x'. See details and examples.

- mask:

  The filename of a mask - function will only be evaluated inside the
  mask.

- maskval:

  An integer specifying the value inside the mask where to apply the
  function. If left blank (the default) then anything above 0.5 will be
  considered inside the mask. This argument only works for mincApply,
  not pMincApply.

- reduce:

  Whether to filter masked voxels from the resulting object

## Value

out The output is a matrix with the same number of rows a the file sizes
and the same number of columns as output by the function that was
applied. Cast into one of the MINC classes to make writing it out with
mincWriteVolume easier.

## Details

mincApply allows one to execute any R function at every voxel of a set
of files. mincApply runs in the current R process; see also pMincApply
and qMincApply for parallel variants. Unless the function to be applied
takes a single argument, which will take on the voxel values at every
voxel, a wrapper function has to be written. The mincApply function also
supports passing, instead of a function, a quoted expression containing
a free variable 'x', again taking on all voxel values; this is
illustrated in the examples. pMincApply does not support this. The
output of this function also has to be a simple vector or scalar. Note
that interpreted R can be very slow. Mindnumbingly slow. It therefore
pays to write optimal functions or, whenever available, use the
optimized equivalents. In other words and to give one example, use
mincLm rather than applying lm, and if lm must be applied, try to use
lm.fit rather than plain lm. Also consider using parallel apply methods.
If you encounter memory issues, it could be due to minc file caching.
Consider trying with the environment variable MINC_FILE_CACHE_MB set to
a small value like 1.

## See also

mincWriteVolume,mincMean,mincSd,mincLm,mincAnova

## Examples

``` r
if (FALSE) { # \dontrun{
testdir <- getRMINCTestData()
gf <- read.csv(file.path(testdir, "test_data_set.csv"))
ma <- mincApply(gf$jacobians_fixed_2, quote(mean(x)))
mincWriteVolume(ma, "means.mnc")

### pMincApply needs a function, not a quoted expression
testFunc <- function(x) c(1, 2)

### run the whole thing in parallel on the local machine
pout <- pMincApply(gf$jacobians_fixed_2, testFunc, local = TRUE, cores = 2)
mincWriteVolume(pout, "pmincOut.mnc")

### run the whole thing as 4 jobs on a batch queue, with a
### batchtools configuration file
pout <- pMincApply(gf$jacobians_fixed_2, testFunc, batches = 4,
                   conf_file = "~/.batchtools.conf.R")
} # }
```
