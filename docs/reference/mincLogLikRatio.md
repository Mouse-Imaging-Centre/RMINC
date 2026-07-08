# run log likelihood ratio tests for different mincLmer objects

Computes the log likelihood ratio of 2 or more voxel-wise lmer calls,
testing the hypothesis that the more complex model better fits the data.
Note that it requires the mixed effects to have been fitted with maximum
likelihood, and not restricted maximum likelihood; in other words, if
you want to use these log likelihood tests, make sure to specify
REML=FALSE in mincLmer.

## Usage

``` r
mincLogLikRatio(...)
```

## Arguments

- ...:

  Two or more mincLmer objects

## Value

the voxel wise log likelihood test. Will have a number of columns
corresponding to the number of inputs -1. Note that it resorts the
inputs from lowest to highest degrees of freedom

## See also

`lmer` and
[`mincLmer`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLmer.md)
for description of lmer and mincLmer.
[`mincFDR`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFDR.md)
for using the False Discovery Rate to correct for multiple comparisons,
and
[`mincWriteVolume`](https://mouse-imaging-centre.github.io/RMINC/reference/mincWriteVolume.md)
for outputting the values to MINC files.

## Examples

``` r
if (FALSE) { # \dontrun{
m1 <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc", REML=F)
m2 <- mincLmer(filenames ~ age + I(age^2) + sex + (age|id),
               data=gf, mask="mask.mnc", REML=F)
m3 <- mincLmer(filenames ~ age + I(age^2) + I(age^3) + sex + (age|id),
               data=gf, mask="mask.mnc", REML=F)
llr <- mincLogLikRatio(m1, m2, m3)
mincFDR(llr)
mincWriteVolume(llr, "m2vsm3.mnc", "m3")
} # }
```
