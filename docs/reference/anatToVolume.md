# Converts a column from an anatLm model into a volume for viewing or saving

Converts a column from an anatLm model into a volume for viewing or
saving

## Usage

``` r
anatToVolume(anat, labelVolume, column, defs = attr(anat, "definitions"))
```

## Arguments

- anat:

  The anatModel

- labelVolume:

  The volume containing the label definitions

- column:

  String indicating which column to turn into a volume

- defs:

  The path to the label definitions

## Value

The volume with the values from the anatModel

## Examples

``` r
if (FALSE) { # \dontrun{
labelVol <- mincArray(mincGetVolume("some-labels.mnc"))
alm <- anatLm(~Sex, gfBasic, vols)
statsvol <- anatToVolume(alm, labelVol, "F.statistic")
mincPlotSliceSeries(anatVol, statsvol, anatLow = 700, anatHigh = 1400,
  low=1, high=10, symmetric = F, begin=50, end=-50)
} # }
```
