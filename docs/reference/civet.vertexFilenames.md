# Assemble vertex files for a CIVET run

Locate the vertex thickness, area, and volume files for a CIVET run

## Usage

``` r
civet.vertexFilenames(
  gf,
  idvar,
  prefix,
  basedir,
  append = TRUE,
  civetVersion = "1.1.9"
)
```

## Arguments

- gf:

  Data Frame with subject information

- idvar:

  column name in gf with subject IDs

- prefix:

  Prefix specified when CIVET was run

- basedir:

  directory where all CIVET output was stored

- append:

  Whether to append the results to the input gf

- civetVersion:

  Version of CIVET

## Value

A data.frame containing left and right thickness, area, and volume
files.
