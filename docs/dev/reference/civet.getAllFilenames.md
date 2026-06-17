# civet.getAllFilenames

Generates list of filenames output by CIVET

## Usage

``` r
civet.getAllFilenames(
  gf,
  idvar,
  prefix,
  basedir,
  append = TRUE,
  civetVersion = "1.1.9",
  cnf = NULL
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

- cnf:

  A list of configuration information used to parse CIVET files produced
  by CBRAIN. Unnecessary for CIVETs \< 2. Should be read from CBRAIN
  yaml config file, but can be input manually as a list containing
  thickness_method (tlink or laplacian), thickness_kernel (in mm), and
  atlas as either AAL or DKT

## Value

gf is returned with CIVET filenames

## Details

Prior to running, read.csv may be called to generate the input argument
gf. The results will be stored under the column name CIVETFILES either
in the input gf (if append = TRUE) or in a new gf. Currently only CIVET
versions 1.1.9 and 1.1.12 are supported.

## See also

civet.readAllCivetFiles

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
} # }
```
