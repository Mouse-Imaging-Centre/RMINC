\name{civet.checkVersion}
\Rdversion{1.1}
\alias{civet.checkVersion}
\title{
Check for a known Civet Version
}
\description{
   As different versions of Civet may reflect changes in filename or directory structure,
   it is important that the civet.* functions be validated for each new version of Civet.
   This function checks whether a given Civet version number has been validated for use by 
   these routines.
}

\usage{
civet.checkVersion(civetVersion)
}

\arguments{
	\item{civetVersion}{A string specifying the Civet version number, e.g., \dQuote{1.1.7}}.
}

\details{
    If the passed Civet version cannot be validated, a warning message is sent to std output.
}

\value{Nothing is returned.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}


\examples{
\dontrun{
civetVersion <- "1.1.8"
civet.checkVersion(civetVersion)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{internal}


