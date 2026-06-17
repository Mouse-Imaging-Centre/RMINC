# Vertex Mixed Effects Models

Perform linear mixed effects model fitting for vertex data. vertexLmer
should be used the same way as a straight lmer call, except that the
left hand side of the equation contains vertex filenames rather than an
actual response variable.

## Usage

``` r
vertexLmer(
  formula,
  data,
  mask = NULL,
  parallel = NULL,
  REML = TRUE,
  column = 1,
  control = lmerControl(),
  start = NULL,
  verbose = 0L,
  safely = FALSE,
  summary_type = "fixef"
)
```

## Arguments

- formula:

  the lmer formula, filenames go on left hand side

- data:

  the data frame, all items in formula should be in here

- mask:

  the mask within which lmer is solved

- parallel:

  how many processors to run on (default=single processor). Specified as
  a two element vector, with the first element corresponding to the type
  of parallelization, and the second to the number of processors to use.
  For local running set the first element to "local" or "snowfall" for
  back-compatibility, anything else will be run with batchtools see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/pMincApply.md).
  Leaving this argument NULL runs sequentially and may take a long time.

- REML:

  whether to use use Restricted Maximum Likelihood or Maximum Likelihood

- column:

  Which column to treat as the input from vertex files.

- control:

  lmer control function

- start:

  lmer start function

- verbose:

  lmer verbosity control

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

## Details

`vertexLmer`, like its relative
[mincLmer](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincLmer.md)
provides an interface to running linear mixed effects models at every
vertex. Unlike standard linear models testing hypotheses in linear mixed
effects models is more difficult, since the denominator degrees of
freedom are more difficult to determine. RMINC provides estimating
degrees of freedom using the
[`vertexLmerEstimateDF`](https://mouse-imaging-centre.github.io/RMINC/dev/reference/vertexLmerEstimateDF.md)
function. For the most likely models - longitudinal models with a
separate intercept or separate intercept and slope per subject - this
approximation is likely correct. Be careful in using this approximation
if using more complicated random effects structures.

## See also

`lmer` for description of lmer and lmer formulas;
[`mincLm`](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincLm.md)
