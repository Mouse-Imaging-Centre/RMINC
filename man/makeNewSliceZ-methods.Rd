\name{makeNewSliceZ-methods}
\Rdversion{1.1}
\docType{methods}

\alias{makeNewSliceZ}
\alias{makeNewSliceZ-methods}
\alias{makeNewSliceZ,MincVolumeIO-method}

\keyword{methods}

\title{Method to Create an Axial Slice from a MincVolumeIO Object}

\description{
    This method retrieves slice 1 from a MincVolumeIO object, 
    and uses it to instantiate an axial MincSlice object. If an initialization
    vector is provided, it will be used to initialize the slice values, else
    the slice will be initialized to zeros.
}

\section{Methods}{
\describe{
\item{mincVolume = "MincVolumeIO"}{Only one method/signature.}
}}

\usage{
\S4method{makeNewSliceZ}{MincVolumeIO}(mincVolume,initVector)
}

\arguments{
  \item{mincVolume}{
    The name of a MincVolumeIO object.}

  \item{initVector}{
    An initialization vector that must contain the same number of elements as are in the slice.}
}


\details{
    This method fills a very particular need.  Specifically, let's say that we have loaded 
    a bunch of slices into a \code{\link{MincSliceIO}} matrix object, and we do a computation at every
    voxel.  This computation will result in a vector having a length reflecting the number
    of voxels in the slice.  Obviously, we would want to store this result somewhere, so this
    method permits us to create an empty Z-slice to act as a container for the result.  Once
    the result is in a \code{\link{MincSlice}} object, it can then be easily moved into a 
    \code{\link{MincVolumeIO}} object for subsequent writing.

    Note that the reason for having a \code{makeNewSliceZ} method and no \code{newSliceX} or 
    \code{newSliceY} methods is because the \code{\link{MincSliceIO}} objects currently
    only store axial slices.  If a need arises for these other methods, they shall be
    written.
}

\value{A MincSlice object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{putSlice}} \code{\linkS4class{MincSlice}} \code{\linkS4class{MincVolumeIO}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}

# create an empty volume and axial slice
vol <- makeNewVolume("myResultsVolume.mnc", likeTemplate="mni305PET")

# load slice 20 from 2 3D volumes into a MincSliceIO object
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
slicesMatrix <- readBySlice(c(v305, v305), 20)

# compute the mean at each voxel
sliceMean <- apply(slicesMatrix,1,mean)   

# create an axial slice and init with the mean() result
zSlice <- makeNewSliceZ(vol, sliceMean)

# display/visualize the slice
print(zSlice)
plot(zSlice)

# move the results slice into the output volume
# ... as slice 20
vol <- putSlice(zSlice, vol, 20)
}


