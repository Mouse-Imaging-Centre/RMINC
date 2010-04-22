\name{MincVolumeIO-class}
\Rdversion{1.1}
\docType{class}
\alias{MincVolumeIO-class}
\alias{MincVolumeIO}

\title{Class "MincVolumeIO"}

\description{
    This class is the primary class used when performing volume granularity IO.  As such,
    this class holds all essential volume meta-data (via a \code{\linkS4class{MincInfo}} object),
    required display settings, and, of course, all of the real data data itself -- stored in
    a 3 dimensional array.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("MincVolumeIO", data, dim, dimnames, ...)},
    however, these objects should never need to be instantiated directly by the user.  Instead, these
    objects are more commonly created as a result to a call to \code{\link{mincIO.readVolume}}.
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"array"}}
    \item{\code{mincInfo}:}{Object of class \code{"MincInfo"}}
    \item{\code{volumeIntensityRange}:}{Object of class \code{"numeric"}}
    \item{\code{frameNumber}:}{Object of class \code{"numeric"}}
    \item{\code{volumeType}:}{Object of class \code{"character"}}
    \item{\code{colorMap}:}{Object of class \code{"character"}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{array}"}, from data part.
Class \code{"\linkS4class{matrix}"}, by class "array", distance 2, with explicit test and coerce.
Class \code{"\linkS4class{structure}"}, by class "array", distance 2.
Class \code{"\linkS4class{vector}"}, by class "array", distance 3, with explicit coerce.
Class \code{"\linkS4class{vector}"}, by class "array", distance 5, with explicit test and coerce.
}

\section{Methods}{
  \describe{
    \item{mincIO.getSliceX}{\code{signature(mincVolume = "MincVolumeIO")}: ... }
    \item{mincIO.getSliceY}{\code{signature(mincVolume = "MincVolumeIO")}: ... }
    \item{mincIO.getSliceZ}{\code{signature(mincVolume = "MincVolumeIO")}: ... }
    \item{plot}{\code{signature(x = "MincVolumeIO")}: ... }
    \item{print}{\code{signature(x = "MincVolumeIO")}: ... }
    \item{mincIO.putSlice}{\code{signature(mincSlice = "MincSlice", mincVolume = "MincVolumeIO")}: ... }
    \item{show}{\code{signature(object = "MincVolumeIO")}: ... }
    \item{mincIO.writeVolume}{\code{signature(object = "MincVolumeIO")}: ... }
	 }
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincInfo}} \code{\linkS4class{MincSlice}} \code{\linkS4class{MincSliceIO}}
}

\examples{
showClass("MincVolumeIO")
}
\keyword{classes}
