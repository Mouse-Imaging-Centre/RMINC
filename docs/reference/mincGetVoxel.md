# Retrieve Voxel Values

Return the intensity of a given voxel in a set of minc files

## Usage

``` r
mincGetVoxel(filenames, v1, v2 = NULL, v3 = NULL)
```

## Arguments

- filenames:

  paths to the minc files

- v1:

  Either a 3-element vector of voxel coordinates or the first

- v2:

  the second voxel coordinate if not NULL

- v3:

  the third voxel coordinate if not NULL

## Value

Returns a `mincVoxel` object containing a vector of intensities and
attributes specify the voxel and world coordinates of the values.
