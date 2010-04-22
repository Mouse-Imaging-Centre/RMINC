\name{mincIO.writeVolume-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.writeVolume-methods}
\alias{mincIO.writeVolume}
\alias{mincIO.writeVolume,MincVolumeIO-method}

\keyword{methods}

\title{Methods to Write a 3D Volume}

\description{
    Given a MincVolumeIO object, write the data and meta-data to a minc volume.}


\section{Methods}{
\describe{
\item{object = "MincVolumeIO"}{Only one method/signature.}
}}


\usage{
\S4method{mincIO.writeVolume}{MincVolumeIO}(object, filename, clobber=FALSE)
}


\arguments{
  \item{object}{An instantiated MincVolumeIO object.}
  \item{filename}{The filename to be used when writing the volume.  If this argument is not set, then filename
                  will be taken from the MincVolumeIO object.}
  \item{clobber}{Logical argument to specify whether the output file should be overwritten.}
}
\details{
    Given a valid MincVolumeIO object, the object is written to the file system.  If the file already
    exists, then either the file will be over-written (clobber=TRUE), or an error message will
    be presented (clobber=FALSE).
}

\value{This function is run for its side-effect, so no value is returned.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.makeNewVolume}} \code{\link{mincIO.readVolume}} \code{\linkS4class{MincVolumeIO}} 
}

\examples{
\dontrun{
# create an empty volume
vol1 <- mincIO.makeNewVolume("myIcbmTemplateVolume.mnc", likeTemplate="icbm152")
mincIO.writeVolume(vol1)                          # write out as "myIcbmTemplateVolume.mnc"
mincIO.writeVolume(vol1, "writeEmptyVol01.mnc")   # write out as "writeEmptyVol01.mnc"

vol2 <- vol1 +1                                   # modify the volume and overwrite
mincIO.writeVolume(vol2, "writeEmptyVol01.mnc", clobber=TRUE)
}}

