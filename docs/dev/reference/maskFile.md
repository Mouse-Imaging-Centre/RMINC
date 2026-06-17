# Get or set a mask file

Mask files control which voxels to perform operations on

## Usage

``` r
maskFile(x, strict = TRUE)

maskFile(x) <- value
```

## Arguments

- x:

  The minc object

- strict:

  Whether or not to throw an error if the mask does not exist

- value:

  The replacement value for the mask attribute

## Value

The mask name if getting, or the object invisibly if setting

## Functions

- `maskFile(x) <- value`: setter
