\name{MincVoxelIO-class}
\Rdversion{1.1}
\docType{class}
\alias{MincVoxelIO-class}
\alias{MincVoxelIO}

\title{Class "MincVoxelIO"}
\description{
This class permits the loading of volume data at voxel granularity.  That is, rather
than reading an entire 3D volume or slice, this class allows the user to define a specific
voxel location that is to be read across all frames of a 4D volume, or across
all volumes of a series of 3D volumes, or both.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MincVoxelIO", data, nrow, ncol, byrow, dimnames, ...)},
    however, users are not expected to instantiate these objects directly.  More conveniently, these objects are usually
    created by calls to \code{\link{mincIO.readByVoxel}}.
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}}
    \item{\code{mincInfo}:}{Object of class \code{"MincInfo"}}
    \item{\code{voxelCoords}:}{Object of class \code{"numeric"}}
    \item{\code{worldCoords}:}{Object of class \code{"numeric"}}
    \item{\code{filenames}:}{Object of class \code{"character"}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "MincVoxelIO")}: ... }
	 }
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.readByVoxel}}
}

\examples{
showClass("MincVoxelIO")
}
\keyword{classes}


