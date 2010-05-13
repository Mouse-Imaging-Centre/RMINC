\name{mincIO.setProperty-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.setProperty-methods}
\alias{mincIO.setProperty}
\alias{mincIO.setProperty,MincInfo,character-method}
\alias{mincIO.setProperty,MincVolumeIO,character-method}
\alias{mincIO.setProperty,MincSlice,character-method}
\alias{mincIO.setProperty,MincSliceIO,character-method}
\alias{mincIO.setProperty,MincVoxelIO,character-method}
\keyword{methods}

\title{Set a specified property from a mincIO object}

\description{
    Basically, these methods provide \dQuote{setter} functions to slot values within all mincIO objects. Unlike the \dQuote{getter} functions which allow access to all slots, these setters only allow specific slots to be modified, in order to prevent the user from creating an inconsistent mincIO object.
}

\usage{
\S4method{mincIO.setProperty}{MincInfo,character}(mincIOobj, propertyId, value)

\S4method{mincIO.setProperty}{MincVolumeIO,character}(mincIOobj, propertyId, value)

\S4method{mincIO.setProperty}{MincSlice,character}(mincIOobj, propertyId, value)

\S4method{mincIO.setProperty}{MincSliceIO,character}(mincIOobj, propertyId, value)

\S4method{mincIO.setProperty}{MincVoxelIO,character}(mincIOobj, propertyId, value)
}


\arguments{
  \item{mincIOobj}{A mincIO object of interest.}
  \item{propertyId}{A string specifying the property whose value is to be modified.  Only specific properties are permitted to be modified.  Use the \code{str()} function to see the object identifiers (as shown in the examples).}
  \item{value}{The new value to which the property is to be set.  The value type (e.g. numeric, character, etc.) must match the slot type.}
}


\value{An updated mincIO object containing the new property value.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.getProperty}}
    \code{\linkS4class{MincInfo}} \code{\linkS4class{MincVolumeIO}} \code{\linkS4class{MincSlice}}
    \code{\linkS4class{MincSliceIO}} \code{\linkS4class{MincVoxelIO}}
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

# show the raw attributes before update
str(vol)

# set selected properties
mincIO.setProperty(vol,"filename", "billy_bob.mnc")
mincIO.setProperty(vol,"volumeIntensityRange", c(0,123))
mincIO.setProperty(vol,"volumeType", "label")
mincIO.setProperty(vol,"colorMap", "AAL")

# show the raw attributes before update
str(vol)
}


