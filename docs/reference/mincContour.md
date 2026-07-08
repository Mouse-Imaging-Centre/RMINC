# Draw contour lines from a MINC volume

Draw contour lines from a MINC volume

## Usage

``` r
mincContour(volume, dimension = 2, slice = NULL, ...)
```

## Arguments

- volume:

  the output of mincArray

- dimension:

  the dimension (from 1 to 3)

- slice:

  the slice number

- ...:

  other parameters to pass on to
  [`contour`](https://rdrr.io/r/graphics/contour.html)

## Examples

``` r
if (FALSE) { # \dontrun{
mincImage(mincArray(anatVol), slice=100, col=gray.colors(255))
mincContour(mincArray(anatVol), slice=100, add=T, col=rainbow(2), levels=c(1000, 1400))
} # }
```
