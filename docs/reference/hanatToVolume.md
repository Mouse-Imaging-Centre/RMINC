# Converts a column from a hierarchical tree into a volume for viewing or saving

Converts a column from a hierarchical tree into a volume for viewing or
saving

## Usage

``` r
hanatToVolume(anatTree, labelVolume, column)
```

## Arguments

- anatTree:

  The anatomical tree

- labelVolume:

  The volume containing the label definitions

- column:

  String indicating which column to turn into a volume

## Value

The volume with the values from the anatomical tree

## Examples

``` r
if (FALSE) { # \dontrun{
labelVol <- mincArray(mincGetVolume("some-labels.mnc"))
hLm <- hanatLm(~Sex, gfBasic, vols)
statsvol <- hanatToVolume(hLm, labelVol, "F.statistic")
mincPlotSliceSeries(anatVol, statsvol, anatLow = 700, anatHigh = 1400,
  low=1, high=10, symmetric = F, begin=50, end=-50)
} # }
```
