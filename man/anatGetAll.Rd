\name{anatGetAll}
\alias{anatGetAll}

\title{Get values given a set of files and an atlas}

\description{Computes volumes, means, sums, and similar values across a
  segmented atlas}
\usage{
anatGetAll(filenames, atlas, method=method, defs=defs, dropLabels=FALSE, side="both")
}

\arguments{
  \item{filenames}{A vector of filenames (strings) which contain the
  information to be extracted at every structure in the atlas.}
  \item{atlas}{A single filename containing the atlas definitions. This
  MINC volume has to be of the same sampling (sizes and dimension order)
  as the filenames specified in the first argument and use a separate
  integer for each atlas label.}
  \item{method}{A string specifying the way information is to be
  computed at every voxel. See the details section for the possible
  options and what they mean.}
  \item{defs}{A string pointing to the filename containing the label
  definitions. Used to map the integers in the atlas to a proper name
  for the structure and contains additional information for laterality
  of each structure.}
  \item{dropLabels}{Whether to return a value for every structure in the
  defs or just for the ones actually contained in each file.}
  \item{side}{Three choices, "right", "left", and "both" (the default)
  which specify what labels to obtain.}
}

\details{anatGetAll needs a set of files along with an atlas and a set
  of atlas definitions. In the end it will produce one value per label
  in the atlas for each of the input files. How that value is computed
  depends on the methods argument:

  \describe{
    \item{jacobians}{Each file contains log jacobians, and the volume
    for each atlas label is computed by multiplying the jacobian with
    the voxel volume at each voxel.}
    \item{labels}{Each file contains integer labels (i.e. same as the
    atlas). The volume is computed by counting the number of voxels with
    each label and multiplying by the voxel volume.}
    \item{means}{Each file contains an arbitrary number and the mean of
    all voxels inside each label is computed.}
    \item{sums}{Each file contains an aribtrary number and the sum of
    all voxels inside each label is computed.}
    \item{text}{Each file is a comma separated values text file and is
      simply read in.}
  }

}

\value{A matrix with nrows equal to the number of labels in the atlas
  and ncols equal to the number of files.}
\seealso{anatLm,anatCombineStructures}

\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")

# get a list of volumes using the Jacobian determinants
# for each of the label values present in the input files
anatGetAll(filename = gf$jacobian_determinant, atlas = "labels_for_input_file.mnc", method="jacobians", defs = "c57_brain_atlas_labels.csv", dropLabels = FALSE, side = "both")

}}
\keyword{}
