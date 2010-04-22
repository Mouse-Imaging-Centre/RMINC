\name{mincIO.plotSlicePretty}
\Rdversion{1.1}
\alias{mincIO.plotSlicePretty}
\title{     Display a 2D Image Plot}

\description{
This function is called internally by the MincSlice::plot function to display a given
2D image plot.}

\usage{
mincIO.plotSlicePretty(slice, xAxisLabel, yAxisLabel, aspectRatio, colorMap)
}
\arguments{
  \item{slice}{The 2D slice matrix extracted from the data part of a MincSlice object.}
  \item{xAxisLabel}{The desired x-axis label.}
  \item{yAxisLabel}{The desired x-axis label.}
  \item{aspectRatio}{The aspect ratio to be used for the display.}
  \item{colorMap}{The color map to be used for the display.}
}

\details{
     This is an internal function, called by the MincSlice::plot function, to serve as the common
     display engine for 2D slice displays -- regardless of slice orientation.}

\value{A \pkg{lattice} object, containing the plot.}


\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{plot}} \code{\linkS4class{MincSlice}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

