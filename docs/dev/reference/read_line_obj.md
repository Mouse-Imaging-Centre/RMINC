# Read a BIC-obj line file

Parse the BIC obj format for when the object contains lines instead of a
mesh.

## Usage

``` r
read_line_obj(line_obj)
```

## Arguments

- line_obj:

  Path to the object file of interest

## Value

`bic_lines` object, which is a list of matrices, each matrix corresponds
to one line in the object. The matrices are 3xN matrices of world
coordinates.
