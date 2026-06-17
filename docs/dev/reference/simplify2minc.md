# Collate Minc

Helper function to collate the results of a
[mincApplyRCPP](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincApplyRCPP.md)
family
([pMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/pMincApply.md),
[mcMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mcMincApply.md),
and
[qMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/qMincApply.md))
function. Internally it calls
[simplify_masked](https://mouse-imaging-centre.github.io/RMINC/dev/reference/simplify_masked.md)
and then coerces the result to the appropriate volumetric class (e.g
`mincMultiDim`)

## Usage

``` r
simplify2minc(result_list)
```

## Arguments

- result_list:

  The mincApply results to collate.

## Value

an object of class `mincSingleDim`, codemincMultiDim, or codemincList
depending on the dimensions of elements of the input list, see
[simplify_masked](https://mouse-imaging-centre.github.io/RMINC/dev/reference/simplify_masked.md)
for details.
