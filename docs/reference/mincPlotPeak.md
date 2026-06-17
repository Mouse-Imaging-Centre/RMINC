# Plotting of peaks

Plots a slice containing a peak. Optionally plots a graph of that peak
alongside.

## Usage

``` r
mincPlotPeak(
  peak,
  anatomy,
  statistics,
  dim = 2,
  crossCol = "green",
  crossSize = 4,
  plotFunction = NULL,
  ...
)
```

## Arguments

- peak:

  a row from
  [`mincFindPeaks`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFindPeaks.md)

- anatomy:

  a mincArray for the underlying anatomy

- statistics:

  a mincArray for the stats volume

- dim:

  the dimension (1:3)

- crossCol:

  the colour for the cross-hair

- crossSize:

  the size (cex) of the cross-hair

- plotFunction:

  a function which will produce a graph

- ...:

  other details to pass on to
  [`mincPlotAnatAndStatsSlice`](https://mouse-imaging-centre.github.io/RMINC/reference/mincPlotAnatAndStatsSlice.md)

## Examples

``` r
if (FALSE) { # \dontrun{
peaks <- mincFindPeaks(-log10(qvs), "Neonatal:time.to.sac", "pos",
posThreshold=1.3, minDistance=1)
p <- function(peak) {
  gfTiming$voxel <- mincGetWorldVoxel(gfTiming$reljacobians02,
                                      peak["x"], peak["y"], peak["z"])
  qplot(time.to.sac, exp(voxel), data=gfTiming, colour=Neonatal,
        geom="boxplot") + theme_classic()
}
mincPlotPeak(peaks[1,], anatVol, -log10(mincArray(qvs, "Neonatal:time.to.sac")),
             anatLow=700, anatHigh=1400, low=1, high=4, col=heat.colors(244),
             crossCol = "blue", crossSize = 3, plotFunction = p)
} # }
```
