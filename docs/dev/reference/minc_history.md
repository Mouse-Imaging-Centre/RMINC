# Minc History

Retrieve or edit the history of a MINC Volume

## Usage

``` r
minc.get.history(filename)

minc.append.history(filename, new_history = NULL)
```

## Arguments

- filename:

  A path to a minc volume

- new_history:

  A new line to be added to the history defaults to "\[timestamp\]\>\>\>
  Written out by RMINC"

## Value

a character vector with one element per line of history

## Functions

- `minc.get.history()`: retrieve

- `minc.append.history()`: append
