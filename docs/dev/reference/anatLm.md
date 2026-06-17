# Calculates statistics and coefficients for linear model of specified anat structure

Calculates statistics and coefficients for linear model of specified
anat structure

## Usage

``` r
anatLm(formula, data, anat, subset = NULL, weights = NULL)
```

## Arguments

- formula:

  a model formula. The RHS of the formula may contain one term with a
  matrix. If so only the + operator may be used, and only two terms may
  appear on the RHS

- data:

  a data.frame containing variables in formula

- anat:

  an array of atlas labels vs subject data

- subset:

  rows to be used, by default all are used

- weights:

  An optional set of weights to use the regression, must be one per
  subject

## Value

Returns an object containing the R-Squared,value,coefficients,F and t
statistcs that can be passed directly into anatFDR. Additionally has the
attributes for model,stat type and degrees of freedom.

## See also

mincLm,anatLm,anatFDR

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
rmincLm= anatLm(~ Sex,gf,gf$lobeThickness)
} # }
```
