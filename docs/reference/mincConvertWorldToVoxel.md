# World to Voxel

Convert coordinates in the volumes world space (a continuous coordinate
space) to a discrete set of voxel coordinates

## Usage

``` r
mincConvertWorldToVoxel(filename, v1, v2, v3, nearest_voxel = TRUE)
```

## Arguments

- filename:

  A path to a minc file

- v1:

  First world coordinate

- v2:

  Second world coordinate

- v3:

  Third world coordinate

- nearest_voxel:

  Logical whether to round the result to the nearest voxel

## Value

a 3-component numeric vector of voxel coordinates
