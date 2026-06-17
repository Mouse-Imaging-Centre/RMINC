# Get Selected Civet Filenames

Returns either one or more Civet filenames, depending on file type.

## Usage

``` r
civet.getFilenameClassify(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameGrayMatterPve(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameWhiteMatterPve(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameCsfPve(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameStxT1(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameCerebrumMask(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameSkullMask(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameGrayMatterSurfaces(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameWhiteMatterSurfaces(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameMidSurfaces(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenamesCorticalThickness(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  smoothing = "20mm",
  fullPath = TRUE
)

civet.getFilenamesCorticalArea(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  smoothing = "40mm",
  fullPath = TRUE
)

civet.getFilenamesCorticalVolume(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  smoothing = "40mm",
  fullPath = TRUE
)

civet.getFilenameMeanSurfaceCurvature(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameLinearTransform(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
)

civet.getFilenameNonlinearTransform(
  scanID,
  baseDir,
  civetVersion = "1.1.9",
  fullPath = TRUE
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

- fullPath:

  A boolean specifying whether the function is to return either a
  fully-qualified path (TRUE) or just the filename without path (FALSE).

- smoothing:

  A character code indicating the smoothing level used in computing
  thickness, area, or volume e.g. "20mm"

## Value

Either a string or a list is returned, depending on the number of
filenames returned. Specifically, a single filename is returned as a
string, whereas multiple filenames are returned as named lists.

## Details

The purpose of this function is to facilitate writing code requiring
manipulation of Civet products. To this purpose, we have written a
number of convenience functions which, given the type of file desired
and a path to the Civet output directory, are able to determine and
return the actual filename(s).

## Functions

- `civet.getFilenameClassify()`: Tissue classification

- `civet.getFilenameGrayMatterPve()`: gray matter pve

- `civet.getFilenameWhiteMatterPve()`: white matter pve

- `civet.getFilenameCsfPve()`: csf pve

- `civet.getFilenameStxT1()`: Standard to T1 transform

- `civet.getFilenameCerebrumMask()`: brain mask

- `civet.getFilenameSkullMask()`: skull mask

- `civet.getFilenameGrayMatterSurfaces()`: gray matter surfaces

- `civet.getFilenameWhiteMatterSurfaces()`: white matter surfaces

- `civet.getFilenameMidSurfaces()`: mid surfaces

- `civet.getFilenamesCorticalThickness()`: cortical thickness

- `civet.getFilenamesCorticalArea()`: cortical area

- `civet.getFilenamesCorticalVolume()`: cortical volume

- `civet.getFilenameMeanSurfaceCurvature()`: surface curvature

- `civet.getFilenameLinearTransform()`: linear transform

- `civet.getFilenameNonlinearTransform()`: non-linear transform

## Author

Jim Nikelski <nikelski@bic.mni.mcgill.ca>

## Examples

``` r

if (FALSE) { # \dontrun{
library(RMINC)

# set Civet root path and scan-identifier
basePath <- "~/tmp/ADNI/civet/pipeOut"
scanID = "0221-M-AD"

# get the name of the aggregate tissue classification volume
# ... and then read it
classifyVolname <- civet.getFilenameClassify(scanID, basePath)
classifyVol <- mincIO.readVolume(classifyVolname)

# get the left and right gray matter surface filenames and then
# ... print the names
gmSurfName <- civet.getFilenameGrayMatterSurfaces(scanID, basePath)
print(gmSurfName$left)
print(gmSurfName$right)

# get the various transformation file filenames
lin.xfmName <- civet.getFilenameLinearTransform(scanID, basePath)
print(lin.xfmName)
nlin.xfmNames <- civet.getFilenameNonlinearTransform(scanID, basePath)
print(nlin.xfmNames$xfm)    # name of the nlin xfm file
print(nlin.xfmNames$grid)    # name of the nlin grid file
} # }
```
