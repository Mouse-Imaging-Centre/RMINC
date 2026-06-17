# Get or set a likeVolume

LikeVolumes control the dimensions of an output minc file.

## Usage

``` r
likeVolume(x, strict = TRUE)

likeVolume(x) <- value
```

## Arguments

- x:

  The minc object

- strict:

  Whether or not to throw an error if the likeVolume does not exist

- value:

  The replacement value for the likeVolume attribute

## Value

The likeVolume name if getting, or the object invisibly if setting

## Functions

- `likeVolume(x) <- value`: setter
