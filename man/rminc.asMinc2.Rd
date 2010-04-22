\name{rminc.asMinc2}
\Rdversion{1.1}
\alias{rminc.asMinc2}
\title{
Convert Volume to Minc2 Format
}
\description{
   As the minc2 API is only capable of reading minc2 volumes, this
   function converts a given volume to minc2 format, as required.
}

\usage{
rminc.asMinc2(filename, keepName = TRUE)
}

\arguments{
  \item{filename}{minc1/minc2 filename to check/convert.}
  \item{keepName}{If \code{TRUE}, keep the original filename, but place the file in the temp directory.\cr
                  If \code{FALSE}, generate a pseudo-random filename, and place file in the working directory.}
}
\details{
    If the passed file is already in minc2 format, nothing is done, and the original filename is returned.  Iff a conversion is required, the conversion is done via a system call to \command{mincconvert} -- so, the minc programs need to be installed and on your path.
}

\value{Either the original filename, or a new filename pointing to the newly converted minc2 volume.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
     \code{\link{rminc.isMinc1}} \code{\link{rminc.isMinc2}} \code{\link{rminc.isMinc}} \code{\link{rminc.isReadable}}
}

\examples{
\dontrun{
myVolume1 <- "happy_little_minc1_volume.mnc"
myVolume2 <- rminc.asMinc2(myVolume1)
}
}


