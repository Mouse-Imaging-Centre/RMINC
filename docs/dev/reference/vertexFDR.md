# Vertex False Discovery Rates

Takes the output of a minc modelling function and computes False
Discovery Rate thresholds.

## Usage

``` r
vertexFDR(buffer, ...)

# S3 method for class 'vertexMultiDim'
vertexFDR(buffer, ...)

# S3 method for class 'vertexLmer'
vertexFDR(buffer, ...)
```

## Arguments

- buffer:

  The results of a vertexLm type run.

- ...:

  additional parameters to
  [mincFDR](https://mouse-imaging-centre.github.io/RMINC/dev/reference/mincFDR.md)
  like method

## Value

A object of type `mincQvals` with the same number of columns as the
input. Each column now contains the qvalues for each vertex. Areas
outside the mask (if a mask was specified) will be represented by a
value of 1. The result also has an attribute called "thresholds" which
contains the 1, 5, 10, 15, and 20 percent false discovery rate
thresholds.

## Methods (by class)

- `vertexFDR(vertexMultiDim)`: vertexMultiDim

- `vertexFDR(vertexLmer)`: vertexLmer
