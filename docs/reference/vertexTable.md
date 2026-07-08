# Create a table of vertex values

Read files containing vertex data into a matrix

## Usage

``` r
vertexTable(filenames, column = 1)
```

## Arguments

- filenames:

  paths to the vertex data files suported extensions: .csv/.csv.gz -
  will assume it's a comma-separated file with a header everything else
  , assume a space separated file, possibly gzipped

- column:

  -specify the column id (name or number) if input files have multiple
  columns

## Value

a matrix where each \`column\` is a matrix of vertex data corresponding
to a single file
