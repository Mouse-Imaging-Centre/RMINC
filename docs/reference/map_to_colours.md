# Generate a vector of colours from a map

Generate a vector of colours from a map

## Usage

``` r
map_to_colours(
  colour_map,
  colour_range = NULL,
  colour_default = "grey",
  symmetric = NULL,
  labels = FALSE,
  palette = heat.colors(255)
)
```

## Arguments

- colour_map:

  either a vector with a label/measure/statistic for every vertex or a
  character file path pointing to a file with such a vector in a rowwise
  format.

- colour_range:

  a two element numeric vector indicating the min and max values of
  allowable labels/measures/statistics to be includedon the surface

- colour_default:

  The colour given to vertices excluded by colour_range

- symmetric:

  Whether to have a positive and negative colour scale (not yet
  implemented)

- labels:

  Whether or not the colour_map is a set of discrete labels

- palette:

  A palette, AKA look-up-table, providing a linear colour scale for the
  colours in `colour_map`

## Value

A vector of colours with types corresponding to your input palette
