# Compute Native to Stereotactic Rescaling Factor

Returns a single float scalar reflecting a global rescaling factor
needed to transform the native image to stereotactic space.

## Usage

``` r
civet.computeNativeToStxRescalingFactor(
  scanID,
  baseDir,
  civetVersion = "1.1.9"
)
```

## Arguments

- scanID:

  A string specifying the unique scan-id (and thus sub-directory) within
  the Civet root output directory.

- baseDir:

  A string specifying the Civet root output directory. This directory
  will, in turn, contain all of the scanIDs.

- civetVersion:

  An optional string specifying the version of Civet used to create the
  output. This is significant since filenames and directory structures
  may change across difference versions of Civet.

## Value

A scalar float reflecting the rescaling factor is returned.

## Details

XFM files contain information to transform a given volume to a model. In
the case of Civet and rescaling, the XFM contains the rescaling factors
(x,y,z) needed to transform the Native volume to the model, which
currently, is usually the symmetrical icbm-152 model.

This functuon serves to compute a global rescaling factor by reading the
individual x,y,z rescales from the linear XFM, and returning the
product.

Interpretation of rescaling factors:

- rescale \> 1.0 The native brain is *expanded* to fit model

- rescale \< 1.0 The native brain is *reduced* to fit model

## Author

Jim Nikelski <nikelski@bic.mni.mcgill.ca>

## Examples

``` r
if (FALSE) { # \dontrun{
library(RMINC)

# set Civet root path and scan-identifier
basePath <- "~/tmp/ADNI/civet/pipeOut"
scanID = "0221-M-AD"

# compute the global rescaling factor
rescale <- civet.computeNativeToStxRescalingFactor(scanID, basePath)
print(rescale)
} # }
```
