# Anatomy Linear Mixed Effects Model

Fit a linear mixed effects model for each structure in the results of
[anatGetAll](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md).

## Usage

``` r
anatLmer(
  formula,
  data,
  anat,
  REML = TRUE,
  control = lmerControl(),
  verbose = FALSE,
  start = NULL,
  parallel = NULL,
  safely = FALSE,
  summary_type = "fixef",
  weights = NULL,
  resources = list(),
  conf_file = getOption("RMINC_BATCH_CONF")
)
```

## Arguments

- formula:

  the lmer formula, filenames go on left hand side

- data:

  the data frame, all items in formula should be in here

- anat:

  a subject by label matrix of anatomical summaries typically produced
  by
  [anatGetAll](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md)

- REML:

  whether to use use Restricted Maximum Likelihood or Maximum Likelihood

- control:

  lmer control function

- verbose:

  lmer verbosity control

- start:

  lmer start function

- parallel:

  how many processors to run on (default=single processor). Specified as
  a two element vector, with the first element corresponding to the type
  of parallelization, and the second to the number of processors to use.
  For local running set the first element to "local" or "snowfall" for
  back-compatibility, anything else will be run with batchtools see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/pMincApply.md).
  Leaving this argument NULL runs sequentially and may take a long time.

- safely:

  whether or not to wrap the per-voxel lmer code in an exception
  catching block (`tryCatch`), when TRUE this will downgrade errors to
  warnings and return NA for the result.

- summary_type:

  Either one of

  - fixef: default and equivalent to older versions of RMINC, returns
    fixed effect coefficients and t-values

  - ranef: returns random effect coefficients and t-values

  - both: both fixed and random effects

  - anova: return the F-statistic for each fixed effect

  or a function to be used to generate the summary

- weights:

  weights to be applied to each observation

- resources:

  A list of resources to use for the jobs, for example
  ` list(nodes = 1, memory = "8G", walltime = "01:00:00") `. See
  `system.file("parallel/pbs_script.tmpl", package = "RMINC")` and
  `system.file("parallel/sge_script.tmpl", package = "RMINC")` for more
  examples

- conf_file:

  A batchtools configuration file

## Details

`anatLmer`, like its relative
[mincLmer](https://mouse-imaging-centre.github.io/RMINC/reference/mincLmer.md)
provides an interface to running linear mixed effects models at every
vertex. Unlike standard linear models testing hypotheses in linear mixed
effects models is more difficult, since the denominator degrees of
freedom are more difficult to determine. RMINC provides estimating
degrees of freedom using the
[`anatLmerEstimateDF`](https://mouse-imaging-centre.github.io/RMINC/reference/anatLmerEstimateDF.md)
function. For the most likely models - longitudinal models with a
separate intercept or separate intercept and slope per subject - this
approximation is likely correct. Be careful in using this approximations
if using more complicated random effects structures.
