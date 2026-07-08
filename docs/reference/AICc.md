# Compute the Corrected AIC for an object

Corrected AIC for finite sample sizes. Generally recommended over AIC as
it asymptotically approaches AIC as n-\>infinity.

## Usage

``` r
AICc(object, ...)
```

## Arguments

- object:

  An object with an AICc method

- ...:

  additional parameters for methods

## Details

For use with
[mincLm](https://mouse-imaging-centre.github.io/RMINC/reference/mincLm.md)
objects AICc produces a vector of AICc values, one per voxel.
