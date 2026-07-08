# converts a vector index to the voxel indices in MINC

RMINC stores volume data as 1d arrays. This function gives the
corresponding voxel coordinates (in the dimension order of the volume)
for an index into the 1d array.

## Usage

``` r
mincVectorToVoxelCoordinates(volumeFileName, vectorCoord)
```

## Arguments

- volumeFileName:

  the filename of the MINC volume

- vectorCoord:

  the integer array index to convert

## Value

a vector of length 3 containing the MINC indices in volume dimension
order

## Examples

``` r
if (FALSE) { # \dontrun{
index <- mincVectorToVoxelCoordinates("filename.mnc", 345322)
voxel <- mincGetVoxel(gf$filenames, index)
} # }
```
