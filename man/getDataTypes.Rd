\name{getDataTypes}
\Rdversion{1.1}
\alias{getDataTypes}
\title{
    Get Minc2 Data Types
}
\description{
    Returns a data.frame containing the minc2 data types.
}

\usage{
getDataTypes()
}

\details{
    Returns a data.frame comprised of a minc2 data type numeric code and 2 string descriptors.  Use
    of this data.frame permits us to get the descriptive string, given the numeric code.
}
\value{A data.frame containing the minc2 data classes.}

\note{Internal use only.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\examples{
dataTypes <- getDataTypes()
print(dataTypes)
}

\seealso{
    \code{\link{getDataClasses}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}


