# Estimates degrees of freedom

see
[anatLmerEstimateDF](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLmerEstimateDF.md)

## Usage

``` r
hanatLmerEstimateDF(buffer, n = 50)
```

## Arguments

- buffer:

  input hierarchical tree

- n:

  number of structures to use for DF estimation

## Value

input tree with df added

## Details

The volumes inside the anatomical hierarchy and the data must be of the
same length and ordering.

## Examples

``` r
if (FALSE) { # \dontrun{
vols <- addVolumesToHierarchy(hdefs, allvols)
hLm <- hanatLmer(~Sex + (1|ID), gf, vols)
hLm <- hanatLmerEstimateDF(hLm)
hLm <- hanatFDR(hLm)
thresholds(hLm)
} # }
```
