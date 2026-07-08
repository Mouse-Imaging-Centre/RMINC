# Check file step sizes

Check if a collection of minc files all have the same step size.

## Usage

``` r
check_step_sizes(filenames, atlas = NULL, strict = FALSE, tolerance = 1e-05)
```

## Arguments

- filenames:

  A character vector of minc file names

- atlas:

  An optional atlas to compare against the files

- strict:

  Should differences be treated as errors or warnings

- tolerance:

  The tolerance for comparing step sizes, 10e-6 by default

## Value

A 3-component vector of step sizes.
