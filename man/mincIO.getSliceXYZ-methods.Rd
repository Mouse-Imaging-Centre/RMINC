\name{mincIO.getSliceXYZ-methods}
\Rdversion{1.1}
\docType{methods}

\alias{mincIO.getSliceXYZ}
\alias{mincIO.getSliceX-methods}
\alias{mincIO.getSliceY-methods}
\alias{mincIO.getSliceZ-methods}
\alias{mincIO.getSliceX}
\alias{mincIO.getSliceY}
\alias{mincIO.getSliceZ}
\alias{mincIO.getSliceX,MincVolumeIO-method}
\alias{mincIO.getSliceY,MincVolumeIO-method}
\alias{mincIO.getSliceZ,MincVolumeIO-method}

\keyword{methods}

\title{Methods to Get a Slice from a MincVolumeIO Object}

\description{
    This method retrieves a specific slice from a MincVolumeIO object, 
    and uses it to instantiate a MincSlice object.
}

\section{Methods}{
\describe{
\item{mincVolume = "MincVolumeIO"}{Only one method/signature.}
}}

\usage{
\S4method{mincIO.getSliceX}{MincVolumeIO}(mincVolume, sliceNo)

\S4method{mincIO.getSliceY}{MincVolumeIO}(mincVolume, sliceNo)

\S4method{mincIO.getSliceZ}{MincVolumeIO}(mincVolume, sliceNo)
}

\arguments{
  \item{mincVolume}{
    The name of a MincVolumeIO object, most likely created via a call to \code{\link{mincIO.readVolume}}.
}
  \item{sliceNo}{
    All slice numbers are 1-relative, ranging from inferior-to-superior (axial),
     left-to-right (sagittal), and posterior-to-anterior (coronal).
}
}


\details{
    Although the data within the MincVolumeIO object can be manipulated directly via indexing
    into the 3 dimensional (xyz) structure, extracting a single slice permits visualization and
    may make modification somewhat more efficient.
}

\value{A MincSlice object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.putSlice}} \code{\linkS4class{MincSlice}} \code{\linkS4class{MincVolumeIO}} \code{\link{mincIO.readVolume}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# read slice 60 from the 305 volume
v305PET <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
vol <- mincIO.readVolume(v305PET)      # read the volume
s060 <- mincIO.getSliceX(vol, 60)      # extract sagittal slice 60
print(s060)                            # print summary info
plot(s060)                             # visualize it
}



