# Read a MINC file

Load a 3-dimensional MINC2 volume and returns it as a 1D array

## Usage

``` r
mincGetVolume(filename)
```

## Arguments

- filename:

  A string corresponding to the location of the MINC2 file to be read.

## Value

Returns a vector of mincSingleDim class

## See also

mincWriteVolume

## Examples

``` r
if (FALSE) { # \dontrun{
testdir <- getRMINCTestData()
testfile <- mincGetVolume(file.path(testdir, "brain_cut_out.mnc"))
} # }
```
