# Plot a BIC obj

Create a basic plot of a BIC obj in the current rgl device, opening a
new one if necessary. If colour_map is supplied, an overlay is added to
the mesh

## Usage

``` r
# S3 method for class 'bic_obj'
plot(
  x,
  colour_map = NULL,
  colour_range = NULL,
  colour_default = "grey",
  symmetric = FALSE,
  palette = heat.colors(255),
  labels = FALSE,
  colour_bar = TRUE,
  add = FALSE,
  par = list(),
  ...
)
```

## Arguments

- x:

  A `bic_obj` probably created by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)

- colour_map:

  A numeric vector equal in length to the number of vertices in the
  `bic_obj` or the path to a text file with one line per vertex with
  colour information.

- colour_range:

  a two element numeric vector indicating the min and max values of
  allowable labels/measures/statistics to be includedon the surface

- colour_default:

  The colour given to vertices excluded by colour_range

- symmetric:

  Whether to have a positive and negative colour scale (not yet
  implemented)

- palette:

  A palette, AKA look-up-table, providing a linear colour scale for the
  colours in `colour_map`

- labels:

  whether the statistic map should be treated as discrete labels.

- colour_bar:

  whether to draw a colour bar

- add:

  whether or not to add this object to the current rgl device (if
  possible) defaults to opening a new device

- par:

  A list of plot parameters to pass to
  [add_colour_bar](https://mouse-imaging-centre.github.io/RMINC/reference/add_colour_bar.md)

- ...:

  additional arguments to
  [create_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/create_mesh.md)
  including but not limited to colour, specular, and add_normals

## Value

invisibly returns the mesh object
