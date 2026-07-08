# Performs ANOVA on each region specified

Performs ANOVA on each region specified

## Usage

``` r
anatAnova(formula, data = NULL, anat = NULL, subset = NULL)
```

## Arguments

- formula:

  a model formula

- data:

  a data.frame containing variables in formula

- anat:

  an array of atlas labels vs subject data

- subset:

  rows to be used, by default all are used

## Value

Returns an array with the F-statistic for each model specified by
formula with the following attributes

- model design matrix

- stat-type type of statistic used

- df degrees of freedom of each statistic

## See also

mincAnova,vertexAnova

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET",
                           TRUE,"1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
rmincAnova = anatAnova(~ Sex,gf,gf$lobeThickness);
} # }
```
