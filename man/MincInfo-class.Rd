\name{MincInfo-class}
\Rdversion{1.1}
\docType{class}
\alias{MincInfo-class}
\alias{MincInfo}

\title{Class "MincInfo"}
\description{
    The MincInfo class is used to hold essential minc volume information, which
    can then be passed to various functions.  As such, this
    class is contained within a number of other higher-level classes.  While
    this class can be created directly by the user via the \code{\link{mincIO.readMincInfo}} method,
    typically, objects of this class would be instantiated by other higher-level methods. 
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("MincInfo", ...)},
    although direct creation of a MincInfo object by the user should be
    avoided.  For internal use only.
}

\section{Slots}{
	 \describe{
    \item{\code{volumeDataClass}:}{Object of class \code{"numeric"}}
    \item{\code{volumeDataType}:}{Object of class \code{"numeric"}}
    \item{\code{spaceType}:}{Object of class \code{"character"}}
    \item{\code{nDimensions}:}{Object of class \code{"numeric"}}
    \item{\code{dimInfo}:}{Object of class \code{"data.frame"}}
    \item{\code{nFrames}:}{Object of class \code{"numeric"}}
    \item{\code{timeWidths}:}{Object of class \code{"numeric"}}
    \item{\code{timeOffsets}:}{Object of class \code{"numeric"}}
    \item{\code{volumeIntensityRange}:}{Object of class \code{"numeric"}}
    \item{\code{filename}:}{Object of class \code{"character"}}
  }

}
\section{Methods}{
  \describe{
    \item{mincIO.readVolume}{\code{signature(object = "MincInfo")}: ... }
    \item{show}{\code{signature(object = "MincInfo")}: ... }
	 }
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincSlice}} \code{\linkS4class{MincSliceIO}} \code{\linkS4class{MincVolumeIO}}
}

\examples{
showClass("MincInfo")
}

\keyword{classes}



