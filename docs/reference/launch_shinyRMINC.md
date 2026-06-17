# launch a shiny based inspector

This uses shiny to inspect the output of mincLm or mincAnova. Produces
various views, the ability to plot individual voxels, and get the output
of the stats model at that voxel.

## Usage

``` r
launch_shinyRMINC(
  statsoutput,
  anatVol,
  volumes = NULL,
  keepBetas = FALSE,
  plotcolumns = NULL,
  modelfunc = NULL,
  singleStatType = NULL,
  fdr = NULL,
  anatLow = 700,
  anatHigh = 1400
)
```

## Arguments

- statsoutput:

  the output of mincLm, mincAnova, or mincLmer. Alternatively a single
  statistic volume. Must be compatible with
  [mincArray](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md)

- anatVol:

  the anatomical volume on which to display the statistics. Can be
  passed a filename, a MINC volume (from mincGetVol), or a mincArray.

- volumes:

  a matrix or data frame of volumes for plotting

- keepBetas:

  whether to include the beta coefficients

- plotcolumns:

  extra data to be used for plotting

- modelfunc:

  optional modelling function

- singleStatType:

  if passing a single statistic volume to `statsoutput` what statistic
  type is it. Stat-types "b", "t", and "tlmer" will get symmetric colour
  scales. All others will get one sided colour scales.

- fdr:

  An optional mincFDR object for choosing thresholds.

- anatLow:

  The lower threshold value for displaying the underlying anatomy

- anatHigh:

  The upper threshold value for displaying the underlying anatomy

## Examples

``` r
if (FALSE) { # \dontrun{
vs <- mincLm(reljacobians02 ~ sex*treatment, subset(gfs, treatment != "None"))
anatVol <- mincArray(mincGetVolume("anatomyfile.mnc"))
launch_shinyRMINC(vs, anatVol, volumes=gfs$vols,
                  plotcolumns=gfs[,c("sex", "Neonatal")], keepBetas=F)
} # }
```
