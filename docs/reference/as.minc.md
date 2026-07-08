# Coerce to RMINC object

Coerce a relatively simple object to an RMINC known object (currently
`mincList`, `mincSingleDim`, `mincMultiDim`)

## Usage

``` r
as.minc(x)
```

## Arguments

- x:

  the object to coerce

## Value

if x is a known minc type return it, if it is a list, attempt to reduce
it toa minc object via
[simplify2minc](https://mouse-imaging-centre.github.io/RMINC/reference/simplify2minc.md),
otherwise check if the object has columns, if so reclass it as
`mincMultiDim` otherwise reclass it as a `mincSingleDim`
