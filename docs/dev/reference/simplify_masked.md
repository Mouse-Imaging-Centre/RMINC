# Simplify Masked Results

Take the results of an RMINC function that produces masked results and
coerce it to the appropriate result. Scalar elements are converted to a
vector, Vector elements are converted to an array, data.frame elements
are row-bound, and more complex objects are converted to a list.

## Usage

``` r
simplify_masked(result_list)
```

## Arguments

- result_list:

  The list of potentially masked results

## Value

Scalar elements are converted to a vector, Vector elements are converted
to an array, data.frame elements are row-bound, and more complex objects
are converted to a list. Result lists should be all of the same type of
unexpected errors may occur.
