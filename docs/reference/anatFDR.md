# Anatomy False Discovery Rates

Anatomy False Discovery Rates

## Usage

``` r
anatFDR(buffer, ...)

# S3 method for class 'anatModel'
anatFDR(buffer, ...)

# S3 method for class 'anatLmer'
anatFDR(buffer, ...)
```

## Arguments

- buffer:

  the result of an
  [anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
  type call

- ...:

  additional parameters to `mincLmer` like method

## Value

A object of type `mincQvals` with the same number of columns as the
input. Each column now contains the qvalues for each structure. Areas
outside the mask (if a mask was specified) will be represented by a
value of 1. The result also has an attribute called "thresholds" which
contains the 1, 5, 10, 15, and 20 percent false discovery rate

## Methods (by class)

- `anatFDR(anatModel)`: anatModel

- `anatFDR(anatLmer)`: anatLmerModel
