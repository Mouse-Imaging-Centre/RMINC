# Create descriptive statistics across a series of vertex files

This function is used to compute the mean, standard deviation, sum, or
variance of every vertex in a set of vertex files.

## Usage

``` r
vertexMean(filenames, column = 1)

vertexSum(filenames, column = 1)

vertexVar(filenames, column = 1)

vertexSd(filenames, column = 1)
```

## Arguments

- filenames:

  Filenames of the vertex volumes across which to create the descriptive
  statistic.

- column:

  Which column to treat as the input from vertex files.

## Value

- out:

  The output will be a single vector containing as many elements as
  there are vertices in the input files.

## Functions

- `vertexMean()`: mean

- `vertexSum()`: sum

- `vertexVar()`: var

- `vertexSd()`: standard deviation

## See also

vertexLm

## Examples

``` r

if (FALSE) { # \dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# compute the mean at every voxel of all files.
means <- vertexMean(gf$filenames)
} # }
```
