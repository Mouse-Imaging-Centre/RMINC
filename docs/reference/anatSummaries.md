# anatSummaries

These functions are used to compute the mean, standard deviation, sum,
or variance of every region in an anat structure.

## Usage

``` r
anatMean(anat)

anatSum(anat)

anatVar(anat)

anatSd(anat)
```

## Arguments

- anat:

  anat structure.

## Value

out: The output will be a single vector containing as many elements as
there are regions in the input variable.

## Functions

- `anatMean()`: mean

- `anatSum()`: sum

- `anatVar()`: variance

- `anatSd()`: standard deviation

## Examples

``` r
if (FALSE) { # \dontrun{
testdir <- getRMINCTestData()
gf = read.csv(file.path(testdir, "CIVET_TEST.csv"))
gf = civet.getAllFilenames(gf,"ID","TEST",file.path(testdir, "CIVET"),"TRUE","1.1.12")
gf = civet.readAllCivetFiles(file.path(testdir, "AAL.csv"),gf)
vm <- anatMean(gf$lobeThickness)
} # }
```
