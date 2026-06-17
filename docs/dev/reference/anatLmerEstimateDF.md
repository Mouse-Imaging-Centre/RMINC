# Estimate Degrees of freedom for an anatLmer model

Estimate the degrees of freedom for an `anatLmerModel` object produced
by
[anatLmer](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLmer.md).
See
[anatLmer](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLmer.md)
for more details on degrees of freedom estimation for linear mixed
effects models

## Usage

``` r
anatLmerEstimateDF(model, n = 50, verbose = FALSE)
```

## Arguments

- model:

  an `anatLmerModel`

- n:

  number of structures to use for DF estimation

- verbose:

  Whether or not to print progress

## Value

the same model, now with degrees of freedom set
