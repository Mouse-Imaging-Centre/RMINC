# Performs ANOVA on each vertex point specified

Performs ANOVA on each vertex point specified

## Usage

``` r
vertexAnova(formula, data, subset = NULL, column = 1)
```

## Arguments

- formula:

  a model formula

- data:

  a data.frame containing variables in formula

- subset:

  rows to be used, by default all are used

- column:

  Which column to treat as the input from vertex files.

## Value

Returns an array with the F-statistic for each model specified by
formula with the following attributes:

- model:

  design matrix

- filenames:

  minc file names input

- dimensions:

  dimensions of the statistics matrix

- dimnames:

  names of the dimensions for the statistic matrix

- stat-type:

  types of statistic used

- df:

  degrees of freedom of each statistic

## See also

mincAnova,anatAnova

## Examples

``` r
if (FALSE) { # \dontrun{
testdir <- getRMINCTestData()
gf = read.csv(file.path(testdir, "CIVET_TEST.csv"))
gf = civet.getAllFilenames(gf,"ID","TEST",file.path(testdir, "CIVET"),"TRUE","1.1.12")
gf = civet.readAllCivetFiles(file.path(testdir, "AAL.csv"),gf)
result = vertexAnova(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf)
} # }
```
