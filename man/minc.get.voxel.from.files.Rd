\name{minc.get.voxel.from.files}
\alias{minc.get.voxel.from.files}
\title{Get one voxel from an array of files}
\description{Given an array of MINC2 filenames, returns the specified
  voxel for each file}
\usage{
voxel <- minc.get.voxel.from.files(filenames, v1, v2, v3)
}
\arguments{
  \item{filenames}{A list of MINC2 filenames.}
  \item{v1}{Integer voxel coordinate of the first dimension in the
    files.}
  \item{v2}{Integer voxel coordinate of the second dimension in the
    files.}
  \item{v3}{Integer voxel coordinate of the third dimension in the
    files.}
}
\details{This function only makes sense if all files in the list have
  the same sampling. Coordinates are given as voxel coordinates. So if
  the volumes have a zyx dimension order, then v1 corresponds to the z
  voxel coordinate, v2 to the y voxel coordinate, and v3 to the x voxel
  coordinate.}
\seealso{minc.get.volume}

  