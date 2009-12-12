\name{readVolumeX}
\Rdversion{1.1}
\alias{readVolumeX}

\title{
    Common Volume Read Function
}
\description{
    This function is called by the various \code{readVolume} methods to perform the
    actual volume read.
}
\usage{
readVolumeX(mincInfo, frameNo, volumeType, colorMap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mincInfo}{See \code{\link{readVolume}}.}
  \item{frameNo}{See \code{\link{readVolume}}.}
  \item{volumeType}{See \code{\link{readVolume}}.}
  \item{colorMap}{See \code{\link{readVolume}}.}
}
\details{
    This is an internal function and should not be called directly by the user.}

\value{An instantiated MincVolumeIO object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{readVolume}} \code{\linkS4class{MincVolumeIO}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}


