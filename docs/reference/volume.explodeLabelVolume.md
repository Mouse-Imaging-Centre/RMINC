# Explode a Label Volume into its Components

Given a label volume, this function splits the volume by label and then
returns a list() containing a mask volume for each of the labels.

## Usage

``` r
volume.explodeLabelVolume(label_vol, labels = NULL, civetLabels = TRUE)
```

## Arguments

- label_vol:

  A string containing the fully-qualified path to the input label
  volume.

- labels:

  Options vector of label names

- civetLabels:

  A logical variable indicating whether the label volume is using the
  Civet convention with regards to naming tissue type (e.g.,
  0=background, 1=csf, etc). If TRUE, the returned list components are
  named using Civet tissue types (bg, csf, gm. wm), else components are
  simply labelled by label number e.g. (“label_0”, “label_2”, etc.).

## Value

A list is returned with each list item holding a mask volume reflecting
a particular label.

## Author

Jim Nikelski <nikelski@bic.mni.mcgill.ca>
