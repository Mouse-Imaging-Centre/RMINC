# civet.CreateBrainViewROI

Create a brainview file of a specific ROI

## Usage

``` r
civet.CreateBrainViewROI(
  atlasFile,
  atlasVertices,
  region,
  civetVersion = "1.1.12"
)
```

## Arguments

- atlasFile:

  File with Atlas label (string-label(int)) mappings

- atlasVertices:

  File with vertex label(int) mappings

- region:

  Region to plot -\> Must match name in atlas file exactly

- civetVersion:

  Version of CIVET

## Details

Create a .txt file that can be overlaid unto the model template (usually
in the CIVET models directory) Currently only CIVET version 1.1.12 is
supported.

## See also

civet.CreateBrainViewFile

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
civet.CreateBrainViewROI("/tmp/rminctestdata/AAL.csv",
                         "/tmp/rminctestdata/AAL_atlas_left.txt",
                         "Left Insula")
q()
#(The .txt file is written in the working directory and can be viewed via the command)
#brain-view2 $CIVET_DIR/models/surf_reg_model_left.obj ./LeftInsula.txt
} # }
```
