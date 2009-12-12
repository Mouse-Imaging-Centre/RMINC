\name{makeNewVolumeX}
\Rdversion{1.1}
\alias{makeNewVolumeX}
\title{    Create a Volume}
\description{
    Use the information passed in a MincInfo object to create a new minc volume object.
}
\usage{
makeNewVolumeX(mincInfo)
}
\arguments{
  \item{mincInfo}{An instantiated and completely initialized MincInfo object.}
}
\details{This internal routine is called by all \code{makeNewVolume} methods to do the actual
     volume object instantiation.}
\value{An instantiated, empty, MincVolumeIO object that has been initialized to zeros.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{makeNewVolume}} \code{\linkS4class{MincVolumeIO}} \code{\linkS4class{MincInfo}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

