\name{MincSlice-class}
\Rdversion{1.1}
\docType{class}
\alias{MincSlice-class}
\alias{MincSlice}

\title{Class "MincSlice"}

\description{
    The purpose of a MincSlice object is to hold a 2 dimensional slice
    of a volume.  Slices can be in any of the 3 orthogonal orientations.
    Overloaded print and plot functions permit display of slice summary
    information, as well as slice visualization.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("MincSlice", data, nrow, ncol, byrow, dimnames, ...)}
    although users are not expected to instantiate these objects directly.
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}}
    \item{\code{mincInfo}:}{Object of class \code{"MincInfo"}}
    \item{\code{sliceNumber}:}{Object of class \code{"numeric"}}
    \item{\code{orientation}:}{Object of class \code{"character"}}
    \item{\code{sliceIntensityRange}:}{Object of class \code{"numeric"}}
    \item{\code{volumeType}:}{Object of class \code{"character"}}
    \item{\code{colorMap}:}{Object of class \code{"character"}}
    \item{\code{aspectRatio}:}{Object of class \code{"numeric"}}
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
    \item{plot}{\code{signature(x = "MincSlice")}: ... }
    \item{print}{\code{signature(x = "MincSlice")}: ... }
    \item{mincIO.putSlice}{\code{signature(mincSlice = "MincSlice", mincVolume = "MincVolumeIO")}: ... }
    \item{show}{\code{signature(object = "MincSlice")}: ... }
	 }
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincInfo}} \code{\linkS4class{MincSliceIO}} \code{\linkS4class{MincVolumeIO}}
}

\examples{
showClass("MincSlice")
}

\keyword{classes}
