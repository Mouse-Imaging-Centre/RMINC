# Get values given a set of files and an atlas

Computes volumes, means, sums, and similar values across a segmented
atlas

## Usage

``` r
anatGetAll_old(
  filenames,
  atlas,
  method = "jacobians",
  defs = getOption("RMINC_LABEL_DEFINITIONS"),
  dropLabels = TRUE,
  side = "both"
)
```

## Arguments

- filenames:

  A vector of filenames (strings) which contain the information to be
  extracted at every structure in the atlas.

- atlas:

  A single filename containing the atlas definitions. This MINC volume
  has to be of the same sampling (sizes and dimension order) as the
  filenames specified in the first argument and use a separate integer
  for each atlas label.

- method:

  A string specifying the way information is to be computed at every
  voxel. See the details section for the possible options and what they
  mean.

- defs:

  A string pointing to the filename containing the label definitions.
  Used to map the integers in the atlas to a proper name for the
  structure and contains additional information for laterality of each
  structure. See
  [voxel_atlas_defs](https://mouse-imaging-centre.github.io/RMINC/reference/voxel_atlas_defs.md)
  for details.

- dropLabels:

  Whether to return a value for every structure in the defs or just for
  the ones actually contained in each file.

- side:

  Three choices, "right", "left", and "both" (the default) which specify
  what labels to obtain.

## Value

A matrix with ncols equal to the number of labels in the atlas and nrows
equal to the number of files.

## Details

anatGetAll needs a set of files along with an atlas and a set of atlas
definitions. In the end it will produce one value per label in the atlas
for each of the input files. How that value is computed depends on the
methods argument:

- jacobians - Each file contains log jacobians, and the volume for each
  atlas label is computed by multiplying the jacobian with the voxel
  volume at each voxel.

- labels - Each file contains integer labels (i.e. same as the atlas).
  The volume is computed by counting the number of voxels with each
  label and multiplying by the voxel volume.

- means - Each file contains an arbitrary number and the mean of all
  voxels inside each label is computed.

- sums - Each file contains an aribtrary number and the sum of all
  voxels inside each label is computed.

- text - Each file is a comma separated values text file and is simply
  read in.

## See also

anatLm,anatCombineStructures

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
volumes <- anatGetAll(filenames=filenames$absolute_jacobian,
                      atlas="/tmp/rminctestdata/test_segmentation.mnc",
                      method="jacobians",
                      defs="/tmp/rminctestdata/test_defs.csv")
} # }
```
