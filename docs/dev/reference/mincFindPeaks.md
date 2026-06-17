# Finds peak values in a MINC file or volume

Inputs a MINC file or RMINC output and identifies peak voxels (given
some constraints) within

## Usage

``` r
mincFindPeaks(
  inputStats,
  column = 1,
  direction = "both",
  minDistance = NA,
  threshold = NA,
  posThreshold = NA,
  negThreshold = NA
)
```

## Arguments

- inputStats:

  Either an RMINC object or a filename from which peaks are to be
  determined.

- column:

  If an RMINC object then the column name or number

- direction:

  "both", "pos" for positive peaks only, or "neg" for negative peaks
  only

- minDistance:

  minimum distance, in mm, between peaks

- threshold:

  threshold above which to consider peaks

- posThreshold:

  positive threshold (overrides threshold if both are given)

- negThreshold:

  negative threshold (overrides threshold if both are given)

## Value

The set of peaks as a matrix. The matrix is of dimensions ntags x 7; the
first three columns correspond to the first three dimensions in the same
coordinates as used by mincArray; the next three columns are the x, y,
and z coordinates in world space, and the 7th column is the value at
that location.

## Details

Provides an interface to find_peaks, a command line program that comes
with the MINC library.
