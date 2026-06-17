# Creates a hierarchical anatomical tree

Makes a hierarchical tree using ontogeny from Allen Brain Institute

## Usage

``` r
makeMICeDefsHierachical(defs, abijson)
```

## Arguments

- defs:

  CSV file describing anatomical labels. See details about format

- abijson:

  The JSON file from the Allen Brain Institute

## Value

a data.tree object containing the full anatomical hierarchy

## Details

Takes the label definitions from a CSV file (same as the rest of the
anat family of functions) and places it into a tree (using data.tree).
The CSV file must thus have an extra column giving the corresponding
name of each structure in the Allen Institute's nomenclature.

Currently the left and right structures are leaves at the end of the
tree, and so are the first to be combined into the bilateral volumes. It
is possible that in the future keeping hemispheres separate will be an
option.

## Examples

``` r
if (FALSE) { # \dontrun{
abijson <- "allen.json"
defs <- "Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv"
hdefs <- makeMICeDefsHierachical(defs, abijson)
} # }
```
