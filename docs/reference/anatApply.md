# Apply function over anat structure

This function is used to compute an arbitrary function of every region
in an anat structure.

## Usage

``` r
anatApply(vols, grouping = NULL, method = mean, ...)
```

## Arguments

- vols:

  anatomy volumes

- grouping:

  a factor grouping with which to perform operations

- method:

  The function which to apply \[default mean\]

- ...:

  Extra arguments to the function passed to method

## Value

out: The output will be a single vector containing as many elements as
there are regions in the input variable by the number of groupings

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
vm <- anatApply(gf$lobeThickness,gf$Primary.Diagnosis)
} # }
```
