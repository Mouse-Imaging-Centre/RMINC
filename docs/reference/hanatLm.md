# linear model applied to every item in an anatomical tree

See
[anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
for more detail

## Usage

``` r
hanatLm(formula, data, anatTree)
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
hLm <- hanatLm(~Sex, gf, vols)
} # }
```
