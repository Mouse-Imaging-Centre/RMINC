# Brain Montage

Plot the left and right hemispheres of brain colourized by two colour
maps with a colour bar in the middle

## Usage

``` r
obj_montage(
  left_obj,
  right_obj,
  left_map,
  right_map,
  output = NULL,
  colour_map,
  colour_range = NULL,
  colour_default = "grey",
  colour_bar = TRUE,
  labels = FALSE,
  palette = heat.colors(255),
  symmetric = FALSE,
  par = list(),
  vertical = TRUE,
  layout = if (vertical) function() rgl::mfrow3d(3, 2, byrow = TRUE) else function()
    rgl::mfrow3d(2, 3, byrow = FALSE),
  ...,
  plot_corners = c(100, 100, 900, 900),
  zoom = 1,
  add_normals = TRUE,
  colour_title = "",
  close_on_output = TRUE
)
```

## Arguments

- left_obj:

  A `bic_obj` probably created by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)
  with the left hemisphere of a subject's brain

- right_obj:

  A `bic_obj` probably created by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)
  with the right hemisphere of a subject's brain

- left_map:

  a colour map to apply to the left hemisphere see
  [colour_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/colour_mesh.md)
  for details

- right_map:

  a colour map to apply to the right hemisphere see
  [colour_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/colour_mesh.md)
  for details

- output:

  Either NULL or a file path to write the snapshot.

- colour_map:

  either a vector with a label/measure/statistic for every vertex or a
  character file path pointing to a file with such a vector in a rowwise
  format.

- colour_range:

  a two element numeric vector indicating the min and max values of
  allowable labels/measures/statistics to be includedon the surface

- colour_default:

  The colour given to vertices excluded by colour_range

- colour_bar:

  Whether or not to draw a colour bar in the figure

- labels:

  Whether or not the colour_map is a set of discrete labels

- palette:

  A palette, AKA look-up-table, providing a linear colour scale for the
  colours in `colour_map`

- symmetric:

  Whether to have a positive and negative colour scale (not yet
  implemented)

- par:

  A list of plot parameters to pass to
  [add_colour_bar](https://mouse-imaging-centre.github.io/RMINC/reference/add_colour_bar.md)

- vertical:

  Whether to use a vertical (default) or horizontal colour bar.

- layout:

  A function to generate an rgl subscene layout, generated with mfrow3d
  or layout3d. The function should take no arguments (a thunk).

- ...:

  additonal parameters to be passed to
  [create_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/create_mesh.md)

- plot_corners:

  The coordinates in pixels for the top left and bottom right corners of
  the the rgl device. `c(lx, ly, rx, ry)`

- zoom:

  A zoom factor to apply to each subplot. This is the inverse of what
  you might expect for consistency with rgl. zoom \> 1 zooms out, zoom
  \< zooms in.

- add_normals:

  Whether or not to add normals to the surface objects, see
  [create_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/create_mesh.md)
  for details

- colour_title:

  legend title for the colour bar if requested

- close_on_output:

  Whether or not to close the output after taking a snapshot, defaults
  to TRUE

## Value

The subscenes invisibly.

## Details

This function is designed to do a simple 6 angle plot for statistic maps
of the left and right hemispheres of subject's brain. It defaults to
leaving the rgl device open so that you can take a snapshot after
tweaking the angles, or changing the colour bar with
[add_colour_bar](https://mouse-imaging-centre.github.io/RMINC/reference/add_colour_bar.md).
Its other mode is to take a snapshot after it has finished adding the 3d
objects to the scene, in this mode, whether or not to keep the window
open can be configured with close_on_output.
