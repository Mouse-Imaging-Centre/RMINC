\name{mincIO.makeNewVolume-methods}
\Rdversion{1.1}
\docType{methods}
\alias{mincIO.makeNewVolume-methods}
\alias{mincIO.makeNewVolume}
\alias{mincIO.makeNewVolume,character,numeric,numeric,numeric-method}
\alias{mincIO.makeNewVolume,character,character-method}

\keyword{methods}

\title{Methods to Create a New, Empty 3D Minc Volume Object}

\description{
    There exist 2 different ways to get an instantiated MincVolumeIO object: (1) create one
    as a result of loading a pre-existing minc volume, and (2) create an empty volume
    from scratch.  This method deals with the second way.
}


\section{Methods}{
\describe{
\item{filename="character", dimLength="numeric", dimSteps="numeric", dimStarts="numeric"}{Create from scratch,
                 by explicitly specifying all of the dimensional attributes.}
\item{filename="character", likeTemplate="character"}{Create the volume by using one of 3 templates.}
\item{filename="character", likeFile="character"}{Create the volume by specifying the full path 
                 and volume name of a volume whose properties the new volume should match.}
}}


\usage{
\S4method{mincIO.makeNewVolume}{character,numeric,numeric,numeric}(filename, dimLengths, dimSteps, dimStarts)

\S4method{mincIO.makeNewVolume}{character,character}(filename, likeTemplate)

\S4method{mincIO.makeNewVolume}{character,character}(filename, likeFile)
}


\arguments{
  \item{filename}{The filename of the minc volume to be created.  The name is only placed into the MincInfo header.
        While no IO is actually performed as this time, the \code{mincIO.writeVolume} method will use this name at file
        creation time.}
  \item{dimLengths}{A vector of 3 integers specifying the number of elements in the x, y, and z dimensions
        (in that order)}
  \item{dimSteps}{A vector of 3 doubles specifying the step size of elements in the x, y, and z dimensions
        (also in that order)}
  \item{dimStarts}{A vector of 3 doubles specifying the start offsets of elements in the x, y, and z dimensions
        (yet again in that order)}
  \item{likeTemplate}{In order to facilitate volume creation, rather than specifying all of the above arguments, 
        one can specify 1 of 3 volume templates whose characteristics are already known to the system.  The 3
        valid template names are: (1) \code{icbm152}, (2) \code{mni305linear}, and (3) \code{mni305PET}.}
  \item{likeFile}{One also has the option of specifying a \dQuote{like} volume.  The properties of the new volume
        shall then be set to match those of the \dQuote{like} volume.}
}


\value{An instantiated, empty, MincVolumeIO object that has been initialized to zeros.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincVolumeIO}} \code{\link{mincIO.readVolume}} \code{\link{mincIO.writeVolume}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# create from scratch
vol <- mincIO.makeNewVolume("mymakeNewVolume.mnc", dimLengths=c(128,128,80), 
                                    dimSteps=c(1.34, 1.72, 1.50), 
                                    dimStarts=c(-85.76, -126.08, -37.50))
# use a template
vol1 <- mincIO.makeNewVolume("myIcbmTemplateVolume.mnc", likeTemplate="icbm152")
vol2 <- mincIO.makeNewVolume("my305LinTemplateVolume.mnc", likeTemplate="mni305linear")
vol3 <- mincIO.makeNewVolume("my305PetTemplateVolume.mnc", likeTemplate="mni305PET")

\dontrun{
# use a "like" volume
vol <- mincIO.makeNewVolume("myLikeVolume", likeFile="myOwnWackySampling.mnc")
}
}
