# Colourize a mesh

Add colour information to your mesh, either from a vertex atlas like AAL
or from a statistic/measurement map like those produced by CIVET

## Usage

``` r
colour_mesh(
  mesh,
  colour_map,
  colour_range = NULL,
  colour_default = "grey",
  symmetric = NULL,
  labels = FALSE,
  palette = heat.colors(255)
)
```

## Arguments

- mesh:

  [mesh3d](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html)
  object ideally produced by
  [create_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/create_mesh.md)

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

an `obj_mesh` object descended from
[mesh3d](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html), with
added colour information and an additional `legend` element to be used
in building a colour bar
