\name{writeVolumeX}
\Rdversion{1.1}
\alias{writeVolumeX}

\title{Common Volume Write Function}

\description{
This function is called by the various \code{readWrite} methods to perform the
actual volume write.}

\usage{
writeVolumeX(mincVolume, filename)
}

\arguments{
  \item{mincVolume}{See \code{\link{writeVolume}}.}
  \item{filename}{See \code{\link{writeVolume}}.}
}

\details{
This is an internal function and should not be called directly by the user.}

\value{Nuttin.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{writeVolume}} \code{\linkS4class{MincVolumeIO}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

