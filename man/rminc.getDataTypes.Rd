\name{rminc.getDataTypes}
\Rdversion{1.1}
\alias{rminc.getDataTypes}
\title{
    Get Minc2 Data Types
}
\description{
    Returns a data.frame containing the minc2 data types.
}

\usage{
rminc.getDataTypes()
}

\details{
    Returns a data.frame comprised of a minc2 data type numeric code and 2 string descriptors.  Use
    of this data.frame permits us to get the descriptive string, given the numeric code.
}
\value{A data.frame containing the minc2 data classes.}

\note{Internal use only.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\examples{
dataTypes <- rminc.getDataTypes()
print(dataTypes)
}

\seealso{
    \code{\link{rminc.getDataClasses}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}


