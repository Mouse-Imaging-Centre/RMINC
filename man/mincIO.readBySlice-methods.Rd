\name{mincIO.readBySlice-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.readBySlice-methods}
\alias{mincIO.readBySlice}
\alias{mincIO.readBySlice,character,numeric-method}

\keyword{methods}

\title{Methods to Read a Specific Slice Across all Frames or Volumes}

\description{
    The purpose of this function is to read a specific slice either across all frames of
    a single 4D volume, or across a collection of 3D volumes.  Note that you can only do
    one or the other of these functions, i.e. read 4D frames in a single volume or read multiple 
    3D volumes. If you really need to read data across a bunch of 4D volumes, take a look
    at the \code{\link{MincVoxelIO}} class.
}

\section{Methods}{
\describe{
\item{filenames = "character", sliceNumber = "numeric"}{Only one method/signature.}
}}


\usage{
\S4method{mincIO.readBySlice}{character,numeric}(filenames, sliceNumber, ..., volumeType, colorMap)
}


\arguments{
  \item{filenames}{Either a single 4D Minc volume name, or a character vector containing multiple 3D volume names.}
  \item{sliceNumber}{The 1-relative axial slice number to be read.  Note that currently we only allow reading axial slices, thus the lack of a \dQuote{slice orientation} argument.}
  \item{\dots}{Other yet to be defined parameters.}
  \item{volumeType}{Type of volume(s) being read.  Current options include \dQuote{anatomical} and \dQuote{functional};
                    if not specified, reasonable defaults will be assigned (4D = \dQuote{functional}, 3D = \dQuote{anatomical}).
                    Note that this optional argument merely effects visualization.}
  \item{colorMap}{Color map to be used for display.  Any valid R colormaps can be specified; if not specified, 
                  reasonable defaults will be assigned (4D = \dQuote{rainbow}, 3D = \dQuote{gray}).}
}

\details{Type of processing is dependent upon the dimensionality of the input volumes.

\describe{
  \item{One 4D volume}{In this case, we shall assume that the user wants to extract all of the
                     frames, for the given slice, for the specified volume.  An object is
                     returned that contains an MxN matrix, where M=slice voxels, and 
                     N=slice over frames.\cr}

  \item{One or more 3D volume(s)}{Here we assume that the user wants to retrieve the specified slice across
                     a range of 3D volumes.  The returned object contains an MxN matrix in which
                     M=slice voxels, and N=slice over volumes.}
}}

\value{An instantiated MincSliceIO object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.getSliceFromSliceArray}} \code{\linkS4class{MincSliceIO}} 
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# read slice 80 across all frames and print summary
v4d <- file.path(volDir, "functional_4D.mnc")
v4d_SA <- mincIO.readBySlice(v4d, 80)                  # read
print(v4d_SA)                                          # ... print summary info

s080f006 <- mincIO.getSliceFromSliceArray(v4d_SA,6)    # get the slice from the 6th frame
print(s080f006)                                        # ... print summary and visualize

# read a slice from a series of 3D volumes
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
s022_SA <- mincIO.readBySlice(c(v305, v305), 22)       # read slice 22 for both volumes
print(s022_SA)                                         # ... print summary
s022v002 <- mincIO.getSliceFromSliceArray(s022_SA,2)   # get slice 22 for the 2nd volume
print(s022v002)                                        # ... print summary and visualize
}




