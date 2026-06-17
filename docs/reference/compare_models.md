# Compare a set of massively-univariate models

For each response (voxel, vertex, or structure) compare a set of linear
model formulations with the criterion of your choice (e.g. AIC, AICc,
BIC).

## Usage

``` r
compare_models(object, ..., metric = AICc)
```

## Arguments

- object:

  A model or a list of models, typically
  [mincLm](https://mouse-imaging-centre.github.io/RMINC/reference/mincLm.md),
  [vertexLm](https://mouse-imaging-centre.github.io/RMINC/reference/vertexLm.md),
  or
  [anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
  results.

- ...:

  additional models

- metric:

  A function to apply to the models that extracts a result for each
  independent sub-model. Typical choices are
  [AIC](https://rdrr.io/r/stats/AIC.html),
  [AICc](https://mouse-imaging-centre.github.io/RMINC/reference/AICc.md),
  and [BIC](https://rdrr.io/r/stats/AIC.html). Please note that metrics
  are considered such that lower is better (in following AIC). To use a
  positive metric create a wrapper function that performs the negation,
  for example, to use the un-modified log-likelihood you could pass
  ` metric = function(minc_model){ -minc_model[,"logLik"]} `

## Value

A sub-model x n models `model_comparison` matrix with the metric of
interest.
