\name{readByVoxelX}
\Rdversion{1.1}
\alias{readByVoxelX}
\title{
Common Read Voxel Function}

\description{
This function is called by the various \code{readByVoxel} methods to perform the
actual voxel read.
}

\usage{
readByVoxelX(filenames, voxelCoords)
}

\arguments{
  \item{filenames}{See \code{\link{readByVoxel}}.}
  \item{voxelCoords}{See \code{\link{readByVoxel}}.}
}
\details{
This is an internal function and should not be called directly by the user.
}

\value{An instantiated MincVoxelIO object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{readByVoxel}} \code{\linkS4class{MincVoxelIO}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

