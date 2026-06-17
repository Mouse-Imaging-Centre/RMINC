# Read a CBRAIN project in R

Simplified file input for CBRAIN projects. Attempts to use as much
information as possible from the directory structure and internal CBRAIN
config files to find and read the correct CIVET results for CIVET 2.0.0
and 2.1.0

## Usage

``` r
civet.readCBRAIN(
  path,
  prefix,
  subjects = NULL,
  atlas = "AAL",
  civetVersion = "2.1.0",
  readFiles = TRUE,
  readQC = TRUE,
  flatten = TRUE,
  QCDir = "QC",
  columnsToKeep = "subject",
  confFile = Sys.glob(paste0(path, "/*/*.yml"))[1]
)
```

## Arguments

- path:

  Path to the civet project

- prefix:

  The prefix used to create the subject names

- subjects:

  A character vector specifying which subjects to read in. If not
  specified, All files within `path` that don't match QC, scripts,
  stats, or analysis will be read in.

- atlas:

  Either AAL (default), DKT, or a path to the specific atlas used

- civetVersion:

  the version of CIVET used, either 2.1.0 (default) or 2.0.0

- readFiles:

  logical whether or not to read the files into R or just generate the
  file names

- readQC:

  logical whether to read and merge the QC results (must be used with
  `flatten`)

- flatten:

  logical whether to convert the CIVET results into a `dplyr` compatible
  data frame or leave it in the legacy format

- QCDir:

  The directory, or vector of directories of where to find QC tables

- columnsToKeep:

  Additional columns from
  [civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md)
  to include in the output

- confFile:

  A configuration file produced by CBRAIN, defaults to the first .yml
  file amongst the subjects.

## Value

A data.frame in the format of
[civet.getAllFilenames](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.getAllFilenames.md)
if `readFiles` is FALSE. A data.frame in the format of
[civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md)
if `readFiles` is TRUE and `flatten` is FALSE. And a data.frame in the
format of
[civet.flattenForDplyr](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.flattenForDplyr.md)
if `readQC`, `readFiles`, and `flatten` are all TRUE (default)

## See also

[civet.getAllFilenames](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.getAllFilenames.md)
[civet.readAllCivetFiles](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.readAllCivetFiles.md)
[civet.flattenForDplyr](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.flattenForDplyr.md)
