# Calculates statistics and coefficients for linear model of specified vertex files

Calculates statistics and coefficients for linear model of specified
vertex files

## Usage

``` r
vertexLm(formula, data, subset = NULL, column = 1)
```

## Arguments

- formula:

  a model formula. The RHS of the formula may contain one term with
  filenames. If so only the + operator may be used, and only two terms
  may appear on the RHS

- data:

  a data.frame containing variables in formula

- subset:

  rows to be used, by default all are used

- column:

  Which column to treat as the input from vertex files.

## Value

Returns an object containing the R-Squared value,beta coefficients, F
and t statistcs that can be passed directly into vertexFDR.

## See also

mincLm,anatLm,vertexFDR

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
result = vertexLm(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf)
} # }
```
