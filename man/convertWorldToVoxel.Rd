\name{convertWorldToVoxel}
\Rdversion{1.1}
\alias{convertWorldToVoxel}

\title{
    Convert World to Voxel Coordinates
}
\description{
    Converts 3 dimensional xyz (world) coordinates into 1-relative ijk (voxel) coordinates.
}
\usage{
    convertWorldToVoxel(filename, worldCoords)
}
\arguments{
  \item{filename}{The name of a minc2 volume to be used to do the conversion to voxel space.}
  \item{worldCoords}{A 3 item vector of world coordinates in xyz order.}
}

\details{
    Conversion requires access to a minc2 volume in order to enable the world to voxel mapping.
}
\value{A 3-item vector of 1-relative voxel coordinates.}

\note{
    TODO: I really should be able to pass a MincInfo object and do the conversion within R code,
    without having to jump into the minc2 API.
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{convertVoxelToWorld}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# set the xyz voxel coordinates
xyzCoords <- c(0, 0, 0)

# create a full pathname to the sample volumes
v305PET <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
print(v305PET)

# convert and show it
ijkCoords <- convertWorldToVoxel(v305PET, xyzCoords)
print(ijkCoords)
}



