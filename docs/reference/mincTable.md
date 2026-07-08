# Read a collection of minc volumes

This reads a collection of \`mincVolumes\` into an optionally
file-backed matrix

## Usage

``` r
mincTable(filenames, mask = NULL, file_backed = FALSE, ...)
```

## Arguments

- filenames:

  A character vector of minc filenames. All minc files must have the
  same shape.

- mask:

  Either a character vector with a path to a mask file, a mincVolume, or
  a numeric/logical vector indicating which voxels to include in the
  table.

- file_backed:

  logical, whether to use a file backed matrix for storing the table.

- ...:

  Additional arguments to \`bigstatr::FBM\`

## Value

If \`file_backed\` is \`FALSE\` return a matrix with ncol equal to the
number of files. The number of rows is either the product of the volume
dimensions, or if a mask is supplied, the number of voxels where the
mask \> 0.5. If \`file_backed\` is \`TRUE\` return an \`FBM\` object
from the \`bigstatr\` package.
