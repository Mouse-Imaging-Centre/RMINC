# Get and Set Minc Specific Attributes

Manage auxillary information contained in RMINC objects. Currently
returns all attributes except dimension information

## Usage

``` r
mincAttributes(minc_object)

setMincAttributes(minc_object, updated_attrs)
```

## Arguments

- minc_object:

  the RMINC object of interest typically a `mincMultiDim` or a related
  class

- updated_attrs:

  A named list containing the new attributes to be replaced

## Value

`mincAttributes` returns an attribute list, `setMincAttributes` returns
the updated object

## Functions

- `setMincAttributes()`: setter
