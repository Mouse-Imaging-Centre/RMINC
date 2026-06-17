# a utility function to compute thresholds

a utility function to compute thresholds

## Usage

``` r
mincFDRThresholdVector(
  pvals,
  qvals,
  thresholdFunc = NULL,
  p.thresholds = c(0.01, 0.05, 0.1, 0.15, 0.2)
)
```

## Arguments

- pvals:

  a vector of pvalues

- qvals:

  a vector of corrected qvalues (such as returend by p.adjust)

- thresholdFunc:

  a function that returns the threshold given a vector of pvalues

- p.thresholds:

  the pvalues at which to compute the threshold

  The function should be the quantile function for the distribution
  being tested. For example, for the chi squared distribution the
  function would be: tfunc \<- function(x) qchisq(max(x), df\[\[i\]\],
  lower.tail=FALSE)
