\name{mincIO.putSlice-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.putSlice-methods}
\alias{mincIO.putSlice}
\alias{mincIO.putSlice,MincSlice,MincVolumeIO-method}

\keyword{methods}

\title{Methods to Transfer a Slice into a MincVolumeIO Object}

\description{
    Move a 2D slice, as held within a MincSlice object, into a MincVolumeIO object.
}

\section{Methods}{
\describe{
\item{mincSlice = "MincSlice", mincVolume = "MincVolumeIO"}{Only one method/signature.}
}}


\usage{
\S4method{mincIO.putSlice}{MincSlice,MincVolumeIO}(mincSlice, mincVolume, sliceNo)
}


\arguments{
  \item{mincSlice}{The name of the MincSlice object.}
  \item{mincVolume}{The name of the MincVolumeIO object, to which the slice is to be written.}
  \item{sliceNo}{A slice number, within the MincVolumeIO object, to which the slice is to be written. This 
                 is not a required argument; by default, the slice number recorded within the MincSlice
                 object is used. Thus, only use this argument if you wish to write the slice back
                 to a different slice from which it was originally retrieved.}
  \item{\dots}{Other yet to be defined parameters.}
}

\details{
     The primary purpose of this method is to permit the re-insertion of newly modified slices into a
     MincVolumeIO volume, prior to writing the volume out.  Now of course, the other major use-case would
     consist of creating a new, empty volume with \code{\link{mincIO.makeNewVolume}}, and then inserting newly created
     slices into this volume.  Note that as it's currently not possible to create an empty slice
     (well, not easily, at least), one could simply create an empty volume and then \code{mincIO.getSliceX/Y/Z} 
     any slice, which would then serve as the template empty slice.  Clearly, in this case, the \code{sliceNo}
     argument would need to be specified in order to write the modified template slices to different positions
     within the volume.}

\value{A newly updated MincVolumeIO object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.getSliceXYZ}} \code{\linkS4class{MincSlice}} 
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# modify a slice and then put it back
v305PET <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
vol <- mincIO.readVolume(v305PET)            # read the 3D volume
sx <- mincIO.getSliceX(vol, 60)              # extract an axial slice
sx[10,20] <- 0                               # modify the slice
updatedVol <- mincIO.putSlice(sx, vol)       # put the slice back into the volume
}


