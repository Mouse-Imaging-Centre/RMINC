# Anatomy and Statistics Slice

Plot a given slice through a MINC volume and superimpose statistics on
the slice.

## Usage

``` r
mincPlotAnatAndStatsSlice(
  anatomy,
  statistics,
  slice = NULL,
  dimension = 2,
  low = min(statistics, na.rm = TRUE),
  high = max(statistics, na.rm = TRUE),
  anatLow = min(anatomy, na.rm = TRUE),
  anatHigh = max(anatomy, na.rm = TRUE),
  symmetric = FALSE,
  col = NULL,
  rcol = NULL,
  legend = NULL,
  acol = gray.colors(255, start = 0),
  legendTextColour = "black"
)
```

## Arguments

- anatomy:

  A minc array of the anatomy volume to plot

- statistics:

  optional statistics or label file to overlay on anatomy slices

- slice:

  the voxel index of the slice of interest

- dimension:

  integer denoting which dimension to slice across

- low:

  the minimum statistic to plot

- high:

  the maximum statistic to plot

- anatLow:

  the minimum anatomy intensity to plot

- anatHigh:

  the maximum antomy intensity to plot

- symmetric:

  whether the statistics are symmetric (such as for t-statistics)

- col:

  colours for statistics

- rcol:

  colours for negative statistics if using a symmetric statistic

- legend:

  an optional string to name the legend, indicating desire for a legend

- acol:

  colours to use for the anatomy

- legendTextColour:

  an optional description of the text colour for the legend (or not)

## Value

invisible NULL
