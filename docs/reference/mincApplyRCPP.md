# Perform Arbitrary calculations on a collection of mincVolumes

An RCPP variant of
[mincApply](https://mouse-imaging-centre.github.io/RMINC/reference/mincApply.md),
the primary advantage being that functions of an arbitrary number of
arguments can be passed to mincApplyRCPP.

## Usage

``` r
mincApplyRCPP(
  filenames,
  fun,
  ...,
  mask = NULL,
  maskval = NULL,
  filter_masked = FALSE,
  slab_sizes = c(1, 1, 1),
  return_indices = FALSE,
  collate = simplify2minc
)
```

## Arguments

- filenames:

  The name of the files to apply over

- fun:

  the function to apply

- ...:

  additional parameters to fun

- mask:

  a numeric mask vector

- maskval:

  An integer specifying the value inside the mask where to apply the
  function. If left blank (the default) then anything above 0.5 will be
  considered inside the mask. This argument only works for mincApply,
  not pMincApply.

- filter_masked:

  Whether or not to remove the masked values from the resultant object

- slab_sizes:

  a three element numeric vector indicating the size in voxels of the
  hyperslab to read for each file. Useful for managing memory use -
  larger slabs are faster but require more memory. Sizes must be an even
  factor of their respective volume dimensions.

- return_indices:

  Whether to return the voxel positions of the results generally for
  internal use only.

- collate:

  A function to (potentially) collapse the result list examples include
  linkunlist and [simplify2array](https://rdrr.io/r/base/lapply.html),
  defaulting to
  [simplify2minc](https://mouse-imaging-centre.github.io/RMINC/reference/simplify2minc.md)
  which creates an object of type `mincMultiDim`, `mincSingleDim`, or
  `mincList` depending on the result structure.\
  If you encounter memory issues, it could be due to minc file caching.
  Consider trying with the environment variable MINC_FILE_CACHE_MB set
  to a small value like 1.

## Value

a list of results subject the the collate function
