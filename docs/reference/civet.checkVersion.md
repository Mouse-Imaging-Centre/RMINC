# Check for a known Civet Version

As different versions of Civet may reflect changes in filename or
directory structure, it is important that the civet.\* functions be
validated for each new version of Civet. This function checks whether a
given Civet version number has been validated for use by these routines.

## Usage

``` r
civet.checkVersion(civetVersion)
```

## Arguments

- civetVersion:

  A string specifying the Civet version number, e.g., “1.1.7”

## Value

Nothing is returned.

## Details

If the passed Civet version cannot be validated, a warning message is
sent to std output.

## Author

Jim Nikelski <nikelski@bic.mni.mcgill.ca>

## Examples

``` r
if (FALSE) { # \dontrun{
civetVersion <- "1.1.8"
civet.checkVersion(civetVersion)
} # }
```
