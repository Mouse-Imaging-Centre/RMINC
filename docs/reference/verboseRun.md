# Run function with/without output silenced

used in test bed

## Usage

``` r
verboseRun(expr, verbose = getOption("verbose", FALSE), env = parent.frame())
```

## Arguments

- expr:

  an expression to run

- verbose:

  whether to permit output, or capture it

- env:

  the environment in which to run the expression, defaults to the
  calling environment
