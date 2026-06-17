# Summarize anatGetAll results by hierarchy

Take the results of
[anatGetAll](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatGetAll.md)
and summarize a grouping label, typically the structures group in the
anatomical hierarchy.

## Usage

``` r
anatSummarize(
  anat,
  summarize_by = "hierarchy",
  defs = getOption("RMINC_LABEL_DEFINITIONS"),
  discard_missing = FALSE
)
```

## Arguments

- anat:

  The results of a call to
  [anatGetAll](https://mouse-imaging-centre.github.io/RMINC/dev/reference/anatGetAll.md)

- summarize_by:

  either a data frame with grouping information or a path to label
  definitions. If a data frame is passed it must have two columns:
  "label" containing the structure name (see `colnames(anat)` to check
  the relevant labels) and "grouping" containing the name of the group
  each structure belongs to.

- defs:

  A text file containing the label definitions if `summarize_by` is a
  string

- discard_missing:

  logical controlling how to handle structures with no ("") group
  information. If TRUE filter these structures, if FALSE give each
  structure a group label to match their structure name.

## Value

A matrix with columns representing groups from the hierarchy and rows
representing individuals with values equal to the sum of the individual
members of each group.
