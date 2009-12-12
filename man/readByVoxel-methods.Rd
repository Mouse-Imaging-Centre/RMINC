\name{readByVoxel-methods}
\Rdversion{1.1}
\docType{methods}
\alias{readByVoxel-methods}
\alias{readByVoxel}
\alias{readByVoxel,character-method}

\keyword{methods}

\title{Methods to Read a Single Voxel over Multiple Frames and Volumes}

\description{
     Given a specific x,y,z location in voxel space, return the value
     at that location across all frames and all volumes.
}

\section{Methods}{
\describe{
\item{filenames = "character"}{Only one method/signature.}
}}

\usage{
readByVoxel(filenames, voxelCoords)
\S4method{readByVoxel}{character}(filenames, voxelCoords)
}

\arguments{
  \item{filenames}{One or more minc2 filenames.  All of these files \emph{should} be in the same space, using the same sampling.}
  \item{voxelCoords}{Voxel coordinate to be sampled.  This xyz coordinate triplet should be 1-relative.}
}

\details{
Here we define that class which does IO at the smallest possible level
of granularity: the voxel.  As is the case with the SliceIO class,
the VoxelIO class tries to make it possible to process a great number
of volumes (3D and 4D), without blowing out our memory.  So, it a nutshell, 
this class reads the real value for a specific 
voxel. Furthermore, the user can specify multiple 4D volumes, thus
allowing us to obtain the value of a specific voxel across all frames,
over all volumes.

This function was inspired by a similar function by Jason Lerch, buried 
in the very bowels of the \pkg{RMINC} package.
}

\value{An instantiated \code{\linkS4class{MincVoxelIO}} object.  The dimensionality of the matrix
       embedded within the data part of the object should be [nFiles x nFrames] in size.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{readByVoxelX}} \code{\linkS4class{MincVoxelIO}} 
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# read voxel value at center of Anterior Commissure across 5 volumes
# into a voxel array.  Start with the 3D case.
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
worldCoords <- c(0,0,0)
voxCoords <- convertWorldToVoxel(v305, worldCoords)
voxelArray <- readByVoxel(c(v305, v305, v305, v305, v305), voxCoords)
print(voxelArray)

# 4D case with 34 frames
v4d <- file.path(volDir, "functional_4D.mnc")
worldCoords <- c(0,0,0)
voxCoords <- convertWorldToVoxel(v4d, worldCoords)
voxelArray <- readByVoxel(c(v4d, v4d, v4d, v4d, v4d), voxCoords)
print(voxelArray)
}


