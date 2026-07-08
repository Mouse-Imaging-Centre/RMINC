# Compute the False Discovery Rate for an anatomical hierarchy

Compute the False Discovery Rate for an anatomical hierarchy

## Usage

``` r
hanatFDR(buffer)
```

## Arguments

- buffer:

  What to compute FDR on

## Value

the anatomical hierarchy with FDR added (see details)

## Details

Uses the Benjamini & Yekutieli (2001) algorithm for computing FDR to
account for dependence (which is obviously present in a hierarchy where
parent values are directly derived from their children).

Stores the values using the same name as the incoming statistic but with
a q prepended. I.e. tvalue.Sexmale becomes qtvalue.Sexmale. Note that
the original, input, is modified with these new additions as well.

## Examples

``` r
if (FALSE) { # \dontrun{
vols <- addVolumesToHierarchy(hdefs, allvols)
hLm <- hanatLm(~Sex, gfBasic, vols)
hLm <- hanatFDR(hLm)
thresholds(hLm)
} # }
```
