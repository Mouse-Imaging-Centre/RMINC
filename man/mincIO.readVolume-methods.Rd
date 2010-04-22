\name{mincIO.readVolume-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.readVolume-methods}
\alias{mincIO.readVolume}
\alias{mincIO.readVolume,MincInfo-method}
\alias{mincIO.readVolume,character-method}

\keyword{methods}

\title{Methods to Read a 3D Volume or a Single 4D Frame}

\description{
Reads a 3D volume or a single 4D frame into a 3D array and returns 
a MincVolumeIO object.}


\section{Methods}{
\describe{
\item{filenames = "character"}{Pass a filename.}
\item{filenames = "MincInfo"}{Pass a single MincInfo object.}
}}


\usage{
\S4method{mincIO.readVolume}{character}(object, frameNo = 0, ..., volumeType, colorMap)

\S4method{mincIO.readVolume}{MincInfo}(object, frameNo = 0, ..., volumeType, colorMap)
}

\arguments{
  \item{object}{The first argument should contain either a MincInfo object, or a character string
                containing the name of a valid Minc volume.}
  \item{frameNo}{A 1-relative frame number.  Only required for 4D volumes.}
  \item{\dots}{Other yet to be defined parameters.}
  \item{volumeType}{Type of volume(s) being read.  Current options include \dQuote{anatomical} and \dQuote{functional};
                    if not specified, reasonable defaults will be assigned (4D = \dQuote{functional}, 3D = \dQuote{anatomical}).
                    Note that this optional argument merely effects visualization.}
  \item{colorMap}{Color map to be used for display.  Any valid R colormaps can be specified; if not specified, 
                  reasonable defaults will be assigned (4D = \dQuote{rainbow}, 3D = \dQuote{gray}).}
}

\details{
    This is the primary function to load an entire 3D volume or frame. Not much to say, really.
}

\value{An instantiated MincVolumeIO object.}

\note{Note that since the MincVolumeIO class inherits from \code{array}, the data portion of the object
      is directly accessible via standard array indexing.  This feature therefore would permit the
      user to both retrieve values from and write new values to the volume object.  Be warned, however, that
      due to R's \dQuote{pass-by-value} semantics, \emph{any} modification of the volume object will
      result in the volume being \emph{copied}.  So, if you want to make your machine cry, loop through
      all of the voxels in a volume and change each voxel's value.  The better way to change a lot
      of values within a volume would involve either: (1) use a vectorized function on the entire volume,
      resulting in only a single volume copy, or (2) do your processing by slice, resulting in only the
      slice being copied for each update.
}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{mincIO.makeNewVolume}} \code{\link{mincIO.writeVolume}} \code{\linkS4class{MincVolumeIO}} 
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# 3D volume
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
vol <- mincIO.readVolume(v305)               # load a 3D volume 
print(vol)                                   # ... print and visualize
\dontrun{plot(vol)}

vol <- mincIO.readVolume(v305, colorMap="rainbow")
\dontrun{plot(vol)}                          # show the colormap effect

vol <- mincIO.readVolume(v305, colorMap="heat.colors") 
\dontrun{plot(vol)}                          # ... yet another colormap

# read a frame from a specific dynamic volume
v4d <- file.path(volDir, "functional_4D.mnc")
vol <- mincIO.readVolume(v4d, 10)            # read frame 10
print(vol)                                   # print and visualize it
\dontrun{plot(vol)}                          #
}


