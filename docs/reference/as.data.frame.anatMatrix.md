# Create Anatomy data.frame

Convert an anatomy frame to data.frame for use with tidy tools

## Usage

``` r
# S3 method for class 'anatMatrix'
as.data.frame(x, ...)
```

## Arguments

- x:

  An `anatMatrix` object produced by
  [anatGetAll](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md)

- ...:

  additional parameters for methods

## Value

A data.frame augmented with addition attributes from anat

## Details

Note that structure names are munged such that each structure name has
non-standard characters replaced by underscores
