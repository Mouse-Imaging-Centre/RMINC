# label peaks with the name of the atlas structure they are in

label peaks with the name of the atlas structure they are in

## Usage

``` r
mincLabelPeaks(peaks, atlas, defs = getOption("RMINC_LABEL_DEFINITIONS"))
```

## Arguments

- peaks:

  the output of
  [`mincFindPeaks`](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincFindPeaks.md)

- atlas:

  the atlas volume, either as a
  [`mincArray`](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincArray.md)
  or as a filename

- defs:

  the atlas definitions, same as used for
  [`anatGetAll`](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatGetAll.md)

## Value

the peaks data frame, with an extra column containing the label

## Examples

``` r
if (FALSE) { # \dontrun{
peaks <- mincFindPeaks(-log10(qvs), "Neonatal:time.to.sac",
                       "pos", posThreshold=1.3, minDistance=1)
peaks <- mincLabelPeaks(peaks,
                         atlasVol,
                         defs="Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv")
} # }
```
