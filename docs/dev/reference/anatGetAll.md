# Faster AnatGet

Computes volumes, means, sums, and similar values across a segmented
atlas

## Usage

``` r
anatGetAll(
  filenames,
  atlas = NULL,
  defs = getOption("RMINC_LABEL_DEFINITIONS"),
  method = c("jacobians", "labels", "sums", "means"),
  side = c("both", "left", "right"),
  parallel = NULL,
  strict = TRUE,
  conf_file = getOption("RMINC_BATCH_CONF")
)
```

## Arguments

- filenames:

  A vector of filenames (strings) which contain the information to be
  extracted at every structure in the atlas.

- atlas:

  One of: NULL if the method is labels, a single atlas file for all
  subjects, or a vector of one atlas file per subject.

- defs:

  A string pointing to the filename containing the label definitions.
  Used to map the integers in the atlas to a proper name for the
  structure and contains additional information for laterality of each
  structure. See
  [voxel_atlas_defs](https://mouse-imaging-centre.github.io/RMINC/dev/reference/voxel_atlas_defs.md)
  for details.

- method:

  A string specifying the way information is to be computed at every
  voxel (default is "jacobians"). See the details section for the
  possible options and what they mean.

- side:

  One of three choices, "right", "left", or "both" (the default) which
  specify what labels to obtain.

- parallel:

  how many processors to run on (default=single processor). Specified as
  a two element vector, with the first element corresponding to the type
  of parallelization, and the second to the number of processors to use.
  For local running set the first element to "local" or "snowfall" for
  back-compatibility, anything else will be run with batchtools see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/pMincApply.md)
  Leaving this argument NULL runs sequentially.

- strict:

  check if any files differ in step sizes

- conf_file:

  A batchtools configuration file defaulting to
  `getOption("RMINC_BATCH_CONF")`

## Value

A matrix with ncols equal to the number of labels in the atlas and nrows
equal to the number of files.

## Details

anatGetAll needs a set of files, an atlas (unless method=="labels"), and
a set of atlas definitions. In the end it will produce one value per
label in the atlas for each of the input files. How that value is
computed depends on the methods argument:

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

If multiple atlases are passed in each subject will have the summary
computed with respect to their atlas.
