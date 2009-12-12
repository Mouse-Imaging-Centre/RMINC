\name{MincSliceIO-class}
\Rdversion{1.1}
\docType{class}
\alias{MincSliceIO-class}
\alias{MincSliceIO}

\title{Class "MincSliceIO"}

\description{
    This class permits the loading of volume data at slice granularity.  That is, rather
    than reading an entire 3D volume, this class allows the user to define a specific
    slice number that is to be read either across all frames of a 4D volume, or across
    all volumes of a series of 3D volumes.}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("MincSliceIO", data, nrow, ncol, byrow, dimnames, ...)}, however, users are not expected to instantiate these objects directly.  More conveniently, these objects are usually created by calls to \link{readBySlice}.}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}}
    \item{\code{mincInfo}:}{Object of class \code{"MincInfo"}}
    \item{\code{sliceNumber}:}{Object of class \code{"numeric"}}
    \item{\code{orientation}:}{Object of class \code{"character"}}
    \item{\code{nSlices}:}{Object of class \code{"numeric"}}
    \item{\code{sliceOrigin}:}{Object of class \code{"character"}}
    \item{\code{filenames}:}{Object of class \code{"character"}}
    \item{\code{volumeType}:}{Object of class \code{"character"}}
    \item{\code{colorMap}:}{Object of class \code{"character"}}
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
    \item{getSliceFromSliceArray}{\code{signature(mincSliceMatrix = "MincSliceIO")}: ... }
    \item{print}{\code{signature(x = "MincSliceIO")}: ... }
    \item{show}{\code{signature(object = "MincSliceIO")}: ... }
	 }
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincInfo}} \code{\linkS4class{MincSlice}} \code{\linkS4class{MincVolumeIO}}
}

\examples{
showClass("MincSliceIO")
}
\keyword{classes}
