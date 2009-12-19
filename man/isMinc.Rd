\name{isMinc}
\Rdversion{1.1}
\alias{isMinc}
\title{
    Check File For Minc Volume Format
}
\description{
    This logical function returns \code{TRUE} if the input volume is structured
    as a minc1 (netcdf) or minc2 (hdf5) volume.  Else, \code{FALSE} is returned.
}
\usage{
isMinc(filename)
}

\arguments{
  \item{filename}{Minc filename.}
}

\details{
    This routine calls the \code{\link{isMinc1}} and \code{\link{isMinc2}} functions to do
    the minc checking.
}
\value{ 
    \code{TRUE} is returned for minc volumes, \code{FALSE} for everything else.
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{isMinc1}} \code{\link{isMinc2}} \code{\link{asMinc2}} \code{\link{isReadable}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# check for minc
myFile <- file.path(volDir, "icbm_avg_152_t1_tal_lin.mnc")
if ( isMinc(myFile) ) {
    myFile <- asMinc2(myFile)               # could be minc1, so call converter
    myFileInfo <- readMincInfo(myFile)      # read the mincInfo for this volume
    print(myFileInfo)                       # ... and print it
}
}

