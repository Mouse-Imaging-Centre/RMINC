# Create a table of vertex measures

Read in the vertex data results of a CIVET run

## Usage

``` r
civet.vertexTable(vertex_files)
```

## Arguments

- vertex_files:

  The data frame of vertex file names produced by
  [civet.vertexFilenames](https://mouse-imaging-centre.github.io/RMINC/reference/civet.vertexFilenames.md)

## Value

A named list of matrices, one for each measure (each column of
`vertex_files` other than `ids`). Each matrix has one row for each
vertex and one column for each subject, named by `ids`.
