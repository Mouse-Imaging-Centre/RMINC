# Run a permutation test on a `mincLm` result

Run a permutation test on a `mincLm` result, computing the most extreme
statistic under exchanged response variables. The randomization
distribution of these extremal statistics is returbed.

## Usage

``` r
mincRandomize(
  x,
  R = 500,
  alternative = c("two.sided", "greater"),
  replace = FALSE,
  parallel = NULL,
  columns = grep("tvalue-", colnames(x)),
  resources = list(),
  conf_file = getOption("RMINC_BATCH_CONF")
)

# S3 method for class 'mincLm'
mincRandomize(
  x,
  R = 500,
  alternative = c("two.sided", "greater"),
  replace = FALSE,
  parallel = NULL,
  columns = grep("tvalue-", colnames(x)),
  resources = list(),
  conf_file = getOption("RMINC_BATCH_CONF")
)
```

## Arguments

- x:

  A `mincLm` object.

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

- columns:

  Which columns to compute extrema for, defaults to columns with
  \`tvalue\` in the name.

- resources:

  A list of resources to use for the jobs, for example
  ` list(nodes = 1, memory = "8G", walltime = "01:00:00") `. See
  `system.file("parallel/pbs_script.tmpl", package = "RMINC")` and
  `system.file("parallel/sge_script.tmpl", package = "RMINC")` for more
  examples

- conf_file:

  A batchtools configuration file defaulting to
  `getOption("RMINC_BATCH_CONF")`

## Value

A list with the original object, the randomization distribution of
extremal statistics and configuration args used for computing the
distributions.

## Methods (by class)

- `mincRandomize(mincLm)`: mincLm
