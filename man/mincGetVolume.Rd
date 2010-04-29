\name{mincGetVolume}
\alias{mincGetVolume}

\title{Read a MINC volume into a 1D array}

\description{Load a 3-dimensional MINC2 volume and returns it as a 1D
  array}

\usage{
   mincGetVolume(filename)
}

\arguments{
  \item{filename}{A string corresponding to the location of the MINC2
    file to be read.}
}

\value{Probably returns a vector of mincSingleDim class.}
\keyword{file}
