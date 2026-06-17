# A utility function to give a MINC object spatial dimensions

Currently the plotting functions (mincPlotAnatAndStatsSlice) need an
object with the correct spatial dimensions assigned. While this might in
the future happen at the time those objects are created, for the moment
this utility function works with older style RMINC objects and extract
what you need.

## Usage

``` r
mincArray(volume, dimIndex = 1, maskVolume = NULL)
```

## Arguments

- volume:

  The input volume (from mincLm, mincGetVolume, etc.)

- dimIndex:

  The index into a multidimensional object

- maskVolume:

  The optional mask from which the data was created

## Value

A matrix with 3 dimensions

## Details

At times data is computed within a mask; for example mincTable matrices
can be extracted from masked data. In that case if the mask volume is
passed as an argument then the data will will be reinserted into those
parts of the mask volume greater than zero.

## Note

R uses Fortran indexing, so dimension assignment is c(dim\[3\],
dim\[2\], dim\[1\]) once dimensions are obtained from any libminc
functions (which use C indexing)

## Examples

``` r
if (FALSE) { # \dontrun{
vol <- mincGetVolume("somefile.mnc")
volWithDims <- mincArray(vol)

vs <- mincLm(jacobians ~ genotype, gf)
tvol <- mincArray(vs, 6)
} # }
```
