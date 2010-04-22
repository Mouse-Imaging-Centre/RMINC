\name{mincIO.getSliceFromSliceArray-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.getSliceFromSliceArray-methods}
\alias{mincIO.getSliceFromSliceArray}
\alias{mincIO.getSliceFromSliceArray,MincSliceIO-method}

\keyword{methods}

\title{Methods to Get a Slice from a MincSliceIO Object}

\description{
    This method retrieves a specific slice from a MincSliceIO object, 
    and uses it to instantiate a MincSlice object.
}

\section{Methods}{
\describe{
\item{mincSliceMatrix = "MincSliceIO"}{Only one method/signature.}
}}

\usage{
\S4method{mincIO.getSliceFromSliceArray}{MincSliceIO}(mincSliceMatrix, sliceIndex)
}

\arguments{
  \item{mincSliceMatrix}{
    The name of a MincSliceIO object, most likely created via a call to \code{\link{mincIO.readBySlice}}.
}
  \item{sliceIndex}{
    The column number within the MincSliceIO object data matrix used to create the Slice object. If
    the MincSliceIO object was initialized from a single 4D volume, then each column will represent
    a specific slice across all frames of the dynamic volume.  Thus, in this case, the \code{sliceIndex} 
    would refer to a specific slice at a specific frame.  Similarly, if the MincSliceIO object was 
    initialized from a list of 3D volumes, then the \code{sliceIndex} would refer to a specific slice 
    within a specific volume.
}
}
\details{
    MincSliceIO objects contain their data within a M x N matrix, in which each
    of the N columns store data from a specific slice, obtained from either (1) all frames of a single 4D volume,
    or (2) a collection of 3D volumes. This method can be used to convert one of these matrix
    columns into a MincSlice object for visualization or manipulation.
}

\value{A MincSlice object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.readBySlice}}  \code{\linkS4class{MincSlice}} \code{\linkS4class{MincSliceIO}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}

# read a slice from a dynamic volume
v4d <- file.path(volDir, "functional_4D.mnc")
s080_SA <- mincIO.readBySlice(v4d, 80)                 # read slice 80 across all frames
s080f006 <- mincIO.getSliceFromSliceArray(s080_SA,6)   # get the slice from the 6th frame
print(s080f006)
plot(s080f006)

# read a slice from a series of 3D volumes
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
s020_SA <- mincIO.readBySlice(c(v305, v305), 20)       # read slice 20 for all volumes
s020v002 <- mincIO.getSliceFromSliceArray(s020_SA,2)   # get the slice for the 2nd volume
print(s020v002)
plot(s020v002)
}




