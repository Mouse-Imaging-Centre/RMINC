# linear mixed effects model applied to every item in an anatomical tree

See
[anatLmer](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLmer.md)
for more detail

## Usage

``` r
hanatLmer(formula, data, anatTree, ...)
```

## Arguments

- formula:

  The model formula

- data:

  The data frame containing the right side of the formula

- anatTree:

  The anatomical tree

- ...:

  Extra options (REML, control of parallel execution, etc.) passed on to
  anatLmer

## Value

The anatomical tree with the model information added

## Details

The volumes inside the anatomical hierarchy and the data must be of the
same length and ordering.

## Examples

``` r
if (FALSE) { # \dontrun{
vols <- addVolumesToHierarchy(hdefs, allvols)
hLm <- hanatLmer(~Sex + (1|ID), gf, vols)
} # }
```
