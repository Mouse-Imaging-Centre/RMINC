# Threshold Free Cluster Enhancement

Perform threshold free cluster enhancement as described in Smith and
Nichols (2008). Cluster-like structures are enhanced to allow a hybrid
cluster/voxel analysis to be performed.

## Usage

``` r
mincTFCE(x, ...)

# S3 method for class 'mincSingleDim'
mincTFCE(
  x,
  d = 0.1,
  E = 0.5,
  H = 2,
  side = c("both", "positive", "negative"),
  output_file = NULL,
  keep = is.null(output_file),
  conf_file = getOption("RMINC_BATCH_CONF"),
  ...
)

# S3 method for class 'matrix'
mincTFCE(
  x,
  d = 0.1,
  E = 0.5,
  H = 2,
  side = c("both", "positive", "negative"),
  like_volume,
  ...
)

# S3 method for class 'mincMultiDim'
mincTFCE(
  x,
  d = 0.1,
  E = 0.5,
  H = 2,
  side = c("both", "positive", "negative"),
  like_volume = likeVolume(x),
  ...
)

# S3 method for class 'mincLm'
mincTFCE(
  x,
  R = 500,
  alternative = c("two.sided", "greater"),
  d = 0.1,
  E = 0.5,
  H = 2,
  side = c("both", "positive", "negative"),
  replace = FALSE,
  parallel = NULL,
  resources = list(),
  conf_file = getOption("RMINC_BATCH_CONF"),
  ...
)
```

## Arguments

- x:

  Either a character vector with a single filename, a `mincSingleDim`
  object, or a `matrix` object, or `mincLm` object.

- ...:

  additional arguments for methods

- d:

  The discretization step-size for approximating the threshold integral
  (default .1)

- E:

  The exponent by which to raise the extent statistic (default .5)

- H:

  The exponent by which to raise the height (default 2)

- side:

  Whether to consider positive and negative statistics or both (default
  both)

- output_file:

  A filename for the enhanced volume.

- keep:

  Whether or not to keep the enhanced volume, defaults to whether or not
  a `output_file` was specified.

- conf_file:

  A batchtools configuration file defaulting to
  `getOption("RMINC_BATCH_CONF")`

- like_volume:

  A path to a like volume specifying the dimensions of the output
  volumes

- R:

  number of randomizations to perform

- alternative:

  Whether to consider a one-sided or two-sided alternative hypothesis.
  Default "two-sided", use "greater" for a one sided test.

- replace:

  Sample with or without replacement for the randomization, defaults to
  FALSE (no replacement)

- parallel:

  A two component vector indicating how to parallelize the computation.
  If the first element is "local" the computation will be run via the
  parallel package, otherwise it will be computed using batchtools, see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/pMincApply.md)
  for details. The element should be numeric indicating the number of
  jobs to split the computation into.

- resources:

  A list of resources to use for the jobs, for example
  ` list(nodes = 1, memory = "8G", walltime = "01:00:00") `. See
  `system.file("parallel/pbs_script.tmpl", package = "RMINC")` and
  `system.file("parallel/sge_script.tmpl", package = "RMINC")` for more
  examples

## Value

The behaviour of `mincTFCE` is to perform cluster free enhancement on a
object, in the single dimensional case, a string denoting a minc file or
a `mincSingleDim` object it returns a `mincSingleDim` object with the
result, optionally saving the file if keep is set to true. In the matrix
case each column is converted to a `mincSingleDim` in accordance with
the `likeVolume`, this is then cluster enhanced and recomposed into a
matrix. In the mincLm case a randomization test is performed with the
t-stats enhanced. The return is

- TFCE: A matrix of the tvalue columns after randomization

- randomization_dist: An RxT matrix where R is the number of
  randomizations and T is the number of t-statistic, elements are the
  largest value obtained by the randomized TFCE

- args: Arguments passed to the internal randomzations and TFCE code

## Methods (by class)

- `mincTFCE(mincSingleDim)`: mincSingleDim

- `mincTFCE(matrix)`: matrix

- `mincTFCE(mincMultiDim)`: mincMultiDim

- `mincTFCE(mincLm)`: mincLm
