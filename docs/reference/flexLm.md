# Flexible Linear Model with Voxel-Varying Covariates

Fit the same linear model at every row of a response matrix, allowing
one or more predictors to vary row-by-row in lock-step with the
response. When the formula does not transform any voxel-varying
covariate (no splines, squares, interactions, etc.), `flexLm` detects
this and uses a fast inner-product update of the model matrix instead of
re-evaluating the formula at every row. Otherwise it falls back to
rebuilding the model matrix per row.

## Usage

``` r
flexLm(formula, data, y, ...)
```

## Arguments

- formula:

  The linear model formula. The left-hand side must be the literal name
  `y`; the per-row response vector is bound to that name internally.

- data:

  A data frame containing the model terms that are constant across rows
  of `y`. Has one row per observation (i.e. per column of `y`).

- y:

  A numeric matrix of response values. Rows correspond to
  voxels/vertices/structures; columns to observations. Row `i` of `y` is
  the response vector for the `i`-th model fit.

- ...:

  Named numeric matrices of per-row (e.g. per-voxel) covariates. Each
  must have the same dimensions as `y`, and each argument name must
  match a term referenced in `formula`. At row `i`, row `i` of each
  matrix is substituted into the corresponding column of `data` before
  the model is fit.

## Value

A numeric matrix with `nrow(y)` rows and `2 * p` columns, where `p` is
the number of columns of the model matrix. The first `p` columns are the
regression coefficients (named after the model-matrix columns); the
remaining `p` columns are the corresponding t-statistics, named
`tstatistic:<term>`.

## Examples

``` r
if (FALSE) { # \dontrun{
n_subjects <- 10
n_voxels   <- 100
df         <- data.frame(age = rnorm(n_subjects))
y          <- matrix(rnorm(n_voxels * n_subjects), nrow = n_voxels)
cov_mat    <- matrix(rnorm(n_voxels * n_subjects), nrow = n_voxels)
out <- flexLm(y ~ age + cov, data = df, y = y, cov = cov_mat)
head(out)
} # }
```
