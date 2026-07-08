# mincFDRMask

Returns either the specified mask, the mask associated with the buffer,
or a vector of ones to be used as a mask.

## Usage

``` r
mincFDRMask(mask = NULL, buffer)
```

## Arguments

- mask:

  a mask file or vector to be passed to mincGetMask if left null, the
  buffer is checked for a mask attribute, if no mask is found, a vector
  of ones is used

- buffer:

  a buffer describing a minc volume

## Value

a numeric mask vector
