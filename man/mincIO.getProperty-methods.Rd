\name{mincIO.getProperty-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.getProperty-methods}
\alias{mincIO.getProperty}
\alias{mincIO.getProperty,MincInfo,character-method}
\alias{mincIO.getProperty,MincVolumeIO,character-method}
\alias{mincIO.getProperty,MincSlice,character-method}
\alias{mincIO.getProperty,MincSliceIO,character-method}
\alias{mincIO.getProperty,MincVoxelIO,character-method}
\keyword{methods}

\title{Get a specified property from a mincIO object}

\description{
    Basically, these methods provide accessor functions to slot values within all mincIO objects.  One is still free to access slots directly using the \dQuote{object@slot} syntax, however, use of accessors are to be preferred since they insulate the user from implementation details.
}

\usage{
\S4method{mincIO.getProperty}{MincInfo,character}(mincIOobj, propertyId)

\S4method{mincIO.getProperty}{MincVolumeIO,character}(mincIOobj, propertyId)

\S4method{mincIO.getProperty}{MincSlice,character}(mincIOobj, propertyId)

\S4method{mincIO.getProperty}{MincSliceIO,character}(mincIOobj, propertyId)

\S4method{mincIO.getProperty}{MincVoxelIO,character}(mincIOobj, propertyId)
}


\arguments{
  \item{mincIOobj}{A mincIO object of interest.}
  \item{propertyId}{A string specifying the property whose value is to be returned.  All values that are contained within the object can be retrieved via the appropriate identifier.  Use the \code{str()} function to see the object identifiers (as shown in the examples).}
}


\value{An entity whose structure is appropriate to the value of the property being retrieved.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincInfo}} \code{\linkS4class{MincVolumeIO}} \code{\linkS4class{MincSlice}} \code{\linkS4class{MincSliceIO}} \code{\linkS4class{MincVoxelIO}}
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

# show the raw attributes
str(vol)

# retrieve selected properties
mincIO.getProperty(vol, "sizes")         # dim sizes; in xyz order
mincIO.getProperty(vol, "starts")        # dim starts
mincIO.getProperty(vol, "nDimensions")   # number of dimensions
mincIO.getProperty(vol, "volumeType")    # type of volume
}



