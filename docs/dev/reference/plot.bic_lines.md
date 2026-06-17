# Plot A bic_lines object

Add lines corresponding to the coordinates in a bic_lines object to a
figure. Only draws the projection of the lines on a single dimension, no
regard is given for whether the lines are near the slice of interest.

## Usage

``` r
# S3 method for class 'bic_lines'
plot(x, dimension = 2, ...)
```

## Arguments

- x:

  an `bic_lines` object

- dimension:

  which axis to display the lines on

- ...:

  additional parameters to pass to
  [segments](https://rdrr.io/r/graphics/segments.html)

## Value

NULL invisibly
