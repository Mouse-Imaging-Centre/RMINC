\name{minc.get.volumes}
\alias{minc.get.volumes}
\title{Get data from a list of volumes}
\description{Returns a 2D array of the full volumes of a list of files}
\usage{
volumes <- minc.get.volumes(filenames)
}
\arguments{
  \item{filenames}{A list of filenames containing MINC2 volumes.}
}
\details{Returns a 2D array of n.voxels by n.volumes size. All volumes
  should have the same sampling.}
\seealso{minc.get.volume}
