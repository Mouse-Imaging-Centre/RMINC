# Add volumes to an anatomical hierarchy

Takes a tree of anatomical hierarchies and adds volumes from
[anatGetAll](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md)

## Usage

``` r
addVolumesToHierarchy(hdefs, volumes)
```

## Arguments

- hdefs:

  The anatomical hierarchy

- volumes:

  The matrix of volumes

## Value

The input anatomical hierarchy with volumes added

## Details

Each node in the tree has two new attributes:

- volumes - The matrix of volumes

- meanVolume - The mean of the volumes of that node

Currently propagates volumes up the tree by summing. Future versions
will add propagating through weighted means or other functions

## Examples

``` r
if (FALSE) { # \dontrun{
abijson <- "allen.json"
defs <- "Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv"
hdefs <- makeMICeDefsHierachical(defs, abijson)
allvols <- anatGetAll(gf$filenames)
hanat <- addVolumesToHierarchy(hdefs, allvols)
} # }
```
