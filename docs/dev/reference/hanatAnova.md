# ANOVA applied to every item in an anatomical tree

See
[anatAnova](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatAnova.md)
for more detail

## Usage

``` r
hanatAnova(formula, data, anatTree)
```

## Arguments

- formula:

  The model formula

- data:

  The data frame containing the right side of the formula

- anatTree:

  The anatomical tree

## Value

The anatomical tree with the model information added

## Details

The volumes inside the anatomical hierarchy and the data must be of the
same length and ordering.

## Examples

``` r
if (FALSE) { # \dontrun{
vols <- addVolumesToHierarchy(hdefs, allvols)
hLm <- hanatAnova(~Sex, gf, vols)
} # }
```
