# Create a brain view file

Creates a text file that can be loaded into brain-view2

## Usage

``` r
civet.CreateBrainViewFile(
  dataFile,
  atlasFile,
  atlasVertices,
  outputFileName,
  civetVersion = "1.1.12"
)
```

## Arguments

- dataFile:

  Either the name of a file with atlas labeling or an R array with atlas
  labeling

- atlasFile:

  Text file with map between atlas labels and numbers

- atlasVertices:

  Text file with map between vertex points and atlas numbers

- outputFileName:

  path to file where output will be saved

- civetVersion:

  Version of CIVET used (Default 1.1.12)

## Details

This function will create a txt file that can be loaded into
brain-view2, in order to visualize the results from CIVET. This function
either accepts a text file or R variable as input. If using an R
variable the rows/columns must be labeled with the Atlas Names. The
names are then matched to numbers as given by the AtlasFile, and then
the numbers are matched to vertices in the AtlasVertexFile.

## Examples

``` r
if (FALSE) { # \dontrun{
gf = read.csv("~/SubjectTable.csv")
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET", TRUE, "1.1.12")
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
civet.CreateBrainViewFile(gf$lobeThickness[1,],
                          "/Atlases/AAL/AAL.csv",
                          "/Atlases/AAL/AAL_atlas_left.txt",
                          "leftLobeThickness.txt")
} # }
```
