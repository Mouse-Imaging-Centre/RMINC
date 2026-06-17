# Plot a slice from a MINC volume

Calls the [`image`](https://rdrr.io/r/graphics/image.html) plotting
function from the base R graphics with some additional data munging to
make it easy to work with MINC slices

## Usage

``` r
mincImage(
  volume,
  dimension = 2,
  slice = NULL,
  low = min(volume, na.rm = TRUE),
  high = max(volume, na.rm = TRUE),
  reverse = FALSE,
  underTransparent = FALSE,
  col = gray.colors(255),
  add = FALSE,
  ...
)
```

## Arguments

- volume:

  a MINC volume as returned by
  [`mincArray`](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md)

- dimension:

  the dimension (1-3) to obtain the slice from

- slice:

  the slice number

- low:

  the low end of the range to plot. If not specified then it will be
  estimated based on the volume's histogram.

- high:

  the high end of the range to plot. If not specified then it will be
  estimated based on the volume's histogram.

- reverse:

  whether to look only at negative numbers.

- underTransparent:

  whether to make anything below the low end of the range transparent.

- col:

  a colour palette AKA look-up-table to colourize the slice intensities

- add:

  whether to add the slice to the current plot device or open a new one
  before plotting

- ...:

  other parameters to pass on to the
  [`image`](https://rdrr.io/r/graphics/image.html) function.

## Details

You can get a fuller tutorial on how to use the visualization tools by
executing the following command:
`file.show(system.file("doc/visualizationTutorial.html", package="RMINC"))`

On certain systems the slices are plotted with a reflected y-axis. To
fix this configure `options(RMINC_flip_image = TRUE)`

## Examples

``` r
if (FALSE) { # \dontrun{
mincImage(mincArray(anatVol), slice=100, col=gray.colors(255))
mincImage(mincArray(vs, 6), slice=100, col=rainbow(255),
          underTransparent = T, low=2, high=6, add=T)
} # }
```
