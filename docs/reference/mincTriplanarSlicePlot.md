# Plot Slice Along Each Axis

Show a slice from each axis of minc volume

## Usage

``` r
mincTriplanarSlicePlot(
  anatomy,
  statistics,
  slice = NULL,
  layoutMatrix = NULL,
  ...
)
```

## Arguments

- anatomy:

  a
  [mincArray](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md)
  object containing the source anatomy

- statistics:

  a
  [mincArray](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md)
  object containing a statistic to overlay

- slice:

  3-component vector indicating which slice along each axis

- layoutMatrix:

  A matrix describing the layout for the plots typically produced by
  [layout](https://rdrr.io/r/graphics/layout.html)

- ...:

  extra parameters to be passed to
  [mincPlotAnatAndStatsSlice](https://mouse-imaging-centre.github.io/RMINC/reference/mincPlotAnatAndStatsSlice.md)

## Value

invisible NULL
