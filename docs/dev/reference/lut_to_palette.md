# A tool that returns a color function/palette from color lookup files

A tool that returns a color function/palette from color lookup files

## Usage

``` r
lut_to_palette(
  lookup_table = system.file("luts/spectral", package = "RMINC"),
  alpha = 1
)
```

## Arguments

- lookup_table:

  Either a path to the lookup table file, or the table itself

- alpha:

  A transparency value between 0 and 1 (inclusive)

## Value

A function that takes in an integer specifying the number of colours
required, and returns a vector of colours interpolated from the input
colour lookup file

## Note

For now, the first column of the lookup table (interval where the colour
occurs) is ignored; i.e. equal spacing intervals are assumed

## Examples

``` r
if (FALSE) { # \dontrun{
spectral.colors <- lut_to_palette(system.file("luts/spectral", package = "RMINC"))
spectral.colors(100)
} # }
```
