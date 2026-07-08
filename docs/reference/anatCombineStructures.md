# Combine left and right volumes

Combines left and right labels from volumes obtained from anatGetAll
call

## Usage

``` r
anatCombineStructures(
  vols,
  method = "jacobians",
  defs = getOption("RMINC_LABEL_DEFINITIONS")
)
```

## Arguments

- vols:

  Matrix output from call to anatGetAll

- method:

  A string specifying the way information was computed at every voxel
  ("jacobians","labels","means","sums")

- defs:

  A string pointing to the filename containing the label definitions.
  Used to map the integers in the atlas to a proper name for the
  structure and contains additional information for laterality of each
  structure.

## Value

A matrix with ncols equal to the number of collapsed labels

## Details

anatCombineStructures collapses left and right volume information into
one measure. If "jacobians","sums",or "labels" is selected then the sum
of the left and right is produced, otherwise the mean is produced.

## See also

anatLm,anatGetAll

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
volumes <- anatGetAll(filenames=filenames$absolute_jacobian,
                      atlas="/tmp/rminctestdata/test_segmentation.mnc",
                      method="jacobians",
                      defs="/tmp/rminctestdata/test_defs.csv")
volumes_combined <-
     anatCombineStructures(vols=volumes,
                           method="jacobians",
                           defs="/tmp/rminctestdata/test_defs.csv")
} # }
```
