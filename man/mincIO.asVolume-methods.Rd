\name{mincIO.asVolume-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.asVolume-methods}
\alias{mincIO.asVolume}
\alias{mincIO.asVolume,array,MincVolumeIO,missing-method}
\alias{mincIO.asVolume,array,missing,character-method}
\keyword{methods}

\title{Convert a 3D-Array into a 3D Minc Volume Object}

\description{
    MincVolumeIO objects are basically comprised of a 3D array, plus a number of properties which describe the volume.
    Now, although one can often treat MincVolumeIO just as one would treat any 3D array, nonetheless, there are some computations
    which strip off the class attributes, returning only the 3D array part. The \code{mincIO.asVolume()} function serves to take a
    3D array and convert it back into a MincVolumeIO object.
}

\usage{
\S4method{mincIO.asVolume}{array,MincVolumeIO,missing}(array3D, likeVolObject, likeTemplate)

\S4method{mincIO.asVolume}{array,missing,character}(array3D, likeVolObject, likeTemplate)
}


\arguments{
  \item{array3D}{A 3D array to be converted into a MincVolumeIO object.}
  \item{likeVolObject}{A pre-existing instantiated MincVolumeIO object to serve as a template.}
  \item{likeTemplate}{A character string containing the name of the template volume. See \code{\link{mincIO.makeNewVolume}} for the names and types of templates currently supported.}
}


\value{An instantiated MincVolumeIO object wrapping the passed 3D array.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincVolumeIO}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# load the ICBM-152 volume
vIcbm <- "icbm_avg_152_t1_tal_lin.mnc"
vIcbm <- file.path(volDir, vIcbm)
vol <- mincIO.readVolume(vIcbm)

# make a mask array (not a volume object)
mask3D <- vol[,,]>100000
class(mask3D)

# apply the mask element-wise to the volume
maskedImage <- vol * mask3D
class(maskedImage)

# recast image data into a MincVolumeIO object
vol2 <-mincIO.asVolume(maskedImage, vol)
class(vol2)

# combining previous 2 lines into 1
vol2 <- mincIO.asVolume(vol * mask3D, vol)

}










