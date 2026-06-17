# MINC Slice Series

Plot a series of slices through a minc volume on a given dimension.
Optionally superimpose statistics, and or include a locator contour to
show where slices are.

## Usage

``` r
mincPlotSliceSeries(
  anatomy,
  statistics = NULL,
  dimension = 2,
  mfrow = c(4, 5),
  low = NULL,
  high = NULL,
  anatLow = NULL,
  anatHigh = NULL,
  col = heat.colors(255),
  anatCol = gray.colors(255, start = 0),
  begin = 1,
  end = (dim(anatomy)[dimension] - 1),
  symmetric = FALSE,
  legend = NULL,
  locator = !is.null(legend),
  plottitle = NULL,
  indicatorLevels = NULL,
  discreteStats = FALSE,
  legendHeight = 0.5
)
```

## Arguments

- anatomy:

  A minc array of the anatomy volume to plot

- statistics:

  optional statistics or label file to overlay on anatomy slices

- dimension:

  integer denoting which dimension to slice across

- mfrow:

  A 2 element vector of the form c(rows, columns) indicating the number
  and position of slices to draw - slices are added by rows

- low:

  the minimum statistic to plot, taken from histogram if not supplied
  and not `discreteStats`, otherwise the minimum statistic

- high:

  the maximum statistic to plot, taken from histogram if not supplied
  and not `discreteStats`, otherwise the maximum statistic

- anatLow:

  the minimum anatomy intensity to plot

- anatHigh:

  the maximum antomy intensity to plot

- col:

  colours for statistics or for the anatomy if statistics are not passed

- anatCol:

  colours for the

- begin:

  the first slice to plot, defaults to 1

- end:

  the last slice to plot, defaults to the last slice

- symmetric:

  whether the statistics are symmetric (such as for t-statistics)

- legend:

  an optional string to name the legend, indicating desire for a legend
  (or not)

- locator:

  whether or not to draw the locator, defaults to whether or not you
  requested a legend

- plottitle:

  the title of the plot if desired

- indicatorLevels:

  numeric vector indicating where to draw slice lines on the locator,
  defaults to every slice

- discreteStats:

  Whether stats are discrete values and should should not have their
  range taken from their histogram if unsupplied.

- legendHeight:

  What vertical fraction of the figure should the legend occupy

## Details

You can get a fuller tutorial on how to use the visualization tools by
executing the following command:
`file.show(system.file("doc/visualizationTutorial.html", package="RMINC"))`

On certain systems the slices are plotted with a reflected y-axis. To
fix this configure `options(RMINC_flip_image = TRUE)`

## Examples

``` r
if (FALSE) { # \dontrun{
mincPlotSliceSeries(mincArray(anatVol),           # the anatomical volume
                    mincArray(vs, "tvalue-SexM"), # pull out one column of the stats
                    anatLow=700, anatHigh=1400,   # set anatomy thresholds
                    low=2.5, high=10,             # set stats thresholds
                    symmetric=T,                  # show separate upper and lower
                    begin=25, end=-25  ,          # remove slices from both sides
                    legend="t-statistics")
} # }
```
