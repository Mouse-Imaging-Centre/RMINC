# Convert an anatModel to a data.frame

Convert an anatModel to a data.frame

## Usage

``` r
# S3 method for class 'anatModel'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  structure_column = NULL
)
```

## Arguments

- x:

  An `anatModel` object, for example from
  [anatLm](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
  or
  [anatAnova](https://mouse-imaging-centre.github.io/RMINC/reference/anatAnova.md)

- row.names, optional:

  Passed to [as.data.frame](https://rdrr.io/r/base/as.data.frame.html)

- ...:

  additional parameters for methods

- structure_column:

  If a string, the structure names (the row names of `x`) are also put
  in a column with this name, as the first column. This keeps the names
  when you use tools that drop row names, for example
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html).
  The default `NULL` adds no column.

## Value

A data.frame with one row for each structure

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- anatLm(~ Sex, gf, vols)
as.data.frame(fits, structure_column = "Structure") %>%
  dplyr::filter(`tvalue-SexM` > 3)
} # }
```
