\name{rminc.isMinc2}
\Rdversion{1.1}
\alias{rminc.isMinc2}
\title{
    Test File For Minc2 Volume Format
}
\description{
    This logical function returns \code{TRUE} if the input volume is structured
    as a minc2 (hdf5) volume.  Else, \code{FALSE} is returned.
}
\usage{
rminc.isMinc2(filename)
}

\arguments{
  \item{filename}{Minc filename.}
}

\details{
    This function uses the \command{file} system command to determine whether \code{filename}
    refers to an hdf5 file. As such, execution of the \command{file} command must return
    the string \dQuote{Hierarchical Data Format} for minc2 files.
}
\value{ 
    \code{TRUE} is returned for minc2 volumes, \code{FALSE} for everything else.
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{rminc.isMinc1}} \code{\link{rminc.isMinc}} \code{\link{rminc.asMinc2}} \code{\link{rminc.isReadable}}
}