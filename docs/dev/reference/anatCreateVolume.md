# Create a Statistic Volume

Convert an anatomy volume into a statistic volume with the results of an
[anatLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLm.md)

## Usage

``` r
anatCreateVolume(anat, filename, column = 1)
```

## Arguments

- anat:

  A minc volume object with a `labels` attribute

- filename:

  the filename for the new minc volume

- column:

  which column of the
  [anatLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatLm.md)
  results to use

## Value

Invisibly returns the created volume
