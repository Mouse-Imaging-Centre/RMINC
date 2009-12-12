\name{isReadable}
\Rdversion{1.1}
\alias{isReadable}
\title{
    Check Read Permissions for Files
}
\description{
    Return \code{TRUE} if all passed files are found to be readable.  If \emph{any} file is
    not readable, \code{FALSE} is returned.
}
\usage{
isReadable(filenames)
}

\arguments{
  \item{filenames}{A character vector consisting of one or more filenames.}
}

\details{
    When passed a character vector containing one or more files, the read permissions
    are checked for all files, and \code{TRUE} is returned if all files are
    found to be readable.  This function works for all file types.
}
\value{ 
    \code{TRUE} is returned if all files are readable, else \code{FALSE}.
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\note{
    This code was basically stolen from Jason Lerch's \pkg{RMINC} package.  May dog have mercy on
    my sole.
}

\seealso{
    \code{\link{isMinc}} \code{\link{isMinc1}} \code{\link{isMinc2}} \code{\link{asMinc2}}
}
\examples{
\dontrun{
# check for readability of a bunch of files
myFiles <- list.files(pattern=".mnc$")     # get all minc files in current directory
if ( !isReadable(myFiles) ) {              # test for readability
    error("One of the files was unreadable.  Very sad, that.") 
}
}}

