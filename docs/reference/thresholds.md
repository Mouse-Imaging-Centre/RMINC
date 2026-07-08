# Get Probability Thresholds

Get Probability Thresholds

## Usage

``` r
thresholds(x, ...)

# S3 method for class 'mincQvals'
thresholds(x, ...)

# S3 method for class 'minc_randomization'
thresholds(x, probs = c(0.01, 0.05, 0.1, 0.2), ...)
```

## Arguments

- x:

  A `mincQvals` object, typically computed with `mincFDR` or a
  `minc*_randomzation` type object. methods

- ...:

  extra arguments for methods

- probs:

  What probabilities to compute thresholds for (only applicable with
  randomization objects)

## Value

A matrix of thresholds, accessible with standard matrix indexing

## Methods (by class)

- `thresholds(mincQvals)`: mincQvals

- `thresholds(minc_randomization)`: minc_randomization
