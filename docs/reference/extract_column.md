# Internal function to read in a table or csv file potentially with many columns and extract just the one we need

Internal function to read in a table or csv file potentially with many
columns and extract just the one we need

## Usage

``` r
extract_column(filename, column = 1)
```

## Arguments

- filename:

  path to the vertex data file

- column:

  -specify the column id (name or number) if input files have multiple
  columns

## Value

a vector of vertex data corresponding to the file
