\name{minc.get.hyperslab}
\alias{minc.get.hyperslab}
\title{Get a hyperslab from a MINC2 file}
\usage{
hyperslab <- minc.get.hyperslab("filename.mnc", c(0,0,0),c(10,10,10))
}
\description{Returns a 1D array by extracting a hyperslab of specified
  starts and counts from a MINC2 volume.}
\arguments{
  \item{filename}{The filename of the MINC2 volume from which to extract
    the hyperslab.}
  \item{start}{A 3-dimensional array of voxel coordinate values to
    specify the start of the hyperslab.}
  \item{count}{A 3-dimensional array of voxel coordinate values to
    specify the count of the hyperslab.}
}
\details{This function allows for the extraction of an arbitrary
  contiguous chunk of data from a MINC2 volume. The coordinates are
  voxel coordinates, given in the volume dimension order.}
\seealso{minc.get.volume}
