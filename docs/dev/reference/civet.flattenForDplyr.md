# Flatten CIVET results for dplyr

Convert the data.frame/Matrix/list fusion object produced by
civet.readAllCivetFiles to a data.frame usable with dplyr etc.

## Usage

``` r
civet.flattenForDplyr(civetResults, columnsToKeep)
```

## Arguments

- civetResults:

  A data.frame produced by
  [civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md)

- columnsToKeep:

  vector of column names or indices of columns from the original frame
  to copy to the normalized table. Columns not produced by
  [civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md)
  and not specified here are dropped. See details

## Value

data.frame containing results of
[civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md),
all sub-data.frames and matrices are expanded, non-standard characters
in column names are replaced with underscores, and a prefix denoting the
origin sub-data is given.

## Details

The columnsToKeep vector needs to include the subject identifier or it
will be dropped. Ideally other columns in the vector should be proper
vector/list columns compatible with dplyr
