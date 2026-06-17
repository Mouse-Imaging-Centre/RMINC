# Estimate the degrees of freedom for parameters in a mincLmer model

There is much uncertainty in how to compute p-values for mixed-effects
statistics, related to the correct calculation of the degrees of freedom
of the model (see here <http://glmm.wikidot.com/faq#df>). mincLmer by
default does not return the degrees of freedom as part of its model,
instead requiring an explicit call to a separate function (such as this
one). The implementation here is the Satterthwaite approximation. This
approximation is computed from the data, to avoid the significant
run-time requirement of computing it separate for every voxel, here it
is only computed on a small number of voxels within the mask and the
median DF returned for every variable.

## Usage

``` r
mincLmerEstimateDF(model)
```

## Arguments

- model:

  the output of mincLmer

## Value

the same mincLmer model, now with degrees of freedom set

## See also

[`mincLmer`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLmer.md)
for mixed effects modelling,
[`mincFDR`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFDR.md)
for multiple comparisons corrections.

## Examples

``` r
if (FALSE) { # \dontrun{
vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc")
vs <- mincLmerEstimateDF(vs)
qvals <- mincFDR(vs, mask=attr(vs, "mask"))
qvals
} # }
```
