# selects a few random indices from a volume

Given a filename, select a few random indices using the uniform
distribution from voxels that have a value of 1 (i.e. from a mask
volume)

## Usage

``` r
mincSelectRandomVoxels(volumeFileName, nvoxels = 50, convert = TRUE, ...)
```

## Arguments

- volumeFileName:

  the filename for a MINC volume

- nvoxels:

  the number of voxels to select

- convert:

  whether to convert to MINC voxel space (default) or keep in index
  space

- ...:

  additional arguments

## Value

A vector of length `nvoxels` containing selected voxel indices or if
convert is true a matrix containing the x-y-z coordinates of the
selected voxels.
