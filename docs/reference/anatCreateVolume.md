# Create a Statistic Volume

Convert an anatomy volume into a statistic volume with the results of an
[anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)

## Usage

``` r
anatCreateVolume(anat, filename = NULL, column = 1)
```

## Arguments

- anat:

  A minc volume object with a `labels` attribute

- filename:

  the filename for the new minc volume. If `NULL` (the default), no file
  is written.

- column:

  which column of the
  [anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
  results to use

## Value

The created volume as a
[mincArray](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md).
It is returned invisibly when a file is written. Voxels with labels that
are not in the label definitions are 0.
