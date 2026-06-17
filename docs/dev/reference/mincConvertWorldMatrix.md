# Convert World to Voxel Coordinates

Convert a 3xN matrix of world coordinates to voxel coordinates with
respect to a given minc file.

## Usage

``` r
mincConvertWorldMatrix(filename, world_matrix, nearest_voxel = TRUE)
```

## Arguments

- filename:

  The minc file indicating a coordinate grid

- world_matrix:

  a 3xN matrix of world coordinates

- nearest_voxel:

  logical whether to round the results to the nearest voxel

## Value

a 3xN matrix of voxel coordinates
