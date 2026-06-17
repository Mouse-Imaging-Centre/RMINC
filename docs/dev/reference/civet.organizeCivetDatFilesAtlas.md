# Organizes CIVET .dat files based on an Atlas

Uses an atlas to associate the measurement results in the CIVET .dat
output files with particular structures

## Usage

``` r
civet.organizeCivetDatFilesAtlas(atlasFile, dataFiles, civetVersion = "1.1.12")
```

## Arguments

- atlasFile:

  Character path to a key to the atlas used when running civet. the key
  should be a comma separated file with a header and the following
  form  
  Column 1: Numeric label Column 3: Corresponding structure

- dataFiles:

  character containing paths to .dat files of interest typically
  generated with
  [civet.getAllFilenames](https://mouse-imaging-centre.github.io/RMINC/dev/reference/civet.getAllFilenames.md)

- civetVersion:

  character code for the version of civet used to generate the data
  files
