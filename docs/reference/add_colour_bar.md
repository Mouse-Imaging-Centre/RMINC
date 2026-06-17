# Add a colour bar for a mesh

Add a colour bar that corresponds to the colours in a given colourized
mesh to the current rgl device

## Usage

``` r
add_colour_bar(
  mesh,
  title = "",
  lpos = 0.97,
  rpos = 0.99,
  bpos = NULL,
  tpos = NULL,
  bpos2 = NULL,
  tpos2 = NULL,
  nudge_title_y = 0.5,
  nudge_title_x = 0.82,
  vertical = TRUE,
  ...
)
```

## Arguments

- mesh:

  A `obj_mesh` object created with
  [colour_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/colour_mesh.md)

- title:

  The label to give the colour bar

- lpos:

  the position for the left edge of the colour bar in fraction of the
  plot area, defaults to .97

- rpos:

  the position for the right edge of the colour bar in fraction of the
  plot area, defaults to .99

- bpos:

  the position for the bottom edge of the colour bar in fraction of the
  plot area. In the symmetric case this is for the negative scale.
  Defaults to .25 for non-symmetric and .05 for symmetric.

- tpos:

  the position for the top edge of the colour bar in fraction of the
  plot area. In the symmetric case this is for the negative scale.
  Defaults to .75 for non-symmetric and .45 for symmetric.

- bpos2:

  the position for the bottom edge of the colour bar in fraction of the
  plot area. Only used in the symmetic case for the positive scale,
  defaults to .55

- tpos2:

  the position for the top edge of the colour bar in fraction of the
  plot area. Only used in the symmetric case for the positive scale,
  defaults to .95

- nudge_title_y:

  Offset text from the bottom of a colour bar by y, units are in
  proportions of the colour bar length

- nudge_title_x:

  Offset text from the bottom of a colour bar by x, units are in
  proportions of the colour bar width

- vertical:

  whether to use a vertical or horizontal layout for the colour bar

- ...:

  extra parameters to pass to
  [`plotrix::color.legend`](https://plotrix.github.io/plotrix/reference/color.legend.html)
  and `text`
