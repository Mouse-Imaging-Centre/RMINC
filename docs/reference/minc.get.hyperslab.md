# Get a hyperslab from a MINC2 file

Returns a 1D array by extracting a hyperslab of specified starts and
counts from a MINC2 volume.

## Usage

``` r
minc.get.hyperslab(filename, start, count)
```

## Arguments

- filename:

  The filename of the MINC2 volume from which to extract the

- start:

  A 3-dimensional array of voxel coordinate values to specify the start
  of the hyperslab.

- count:

  A 3-dimensional array of voxel coordinate values to specify the count
  of the hyperslab.

## Value

a numeric vector of size `prod(count)` containing the hyperslab

## Details

This function allows for the extraction of an arbitrary contiguous chunk
of data from a MINC2 volume. The coordinates are voxel coordinates,
given in the volume dimension order.
