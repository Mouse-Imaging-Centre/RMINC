\name{print-methods}
\Rdversion{1.1}
\docType{methods}
\alias{print-methods}
\alias{print,MincInfo-method}
\alias{print,MincSlice-method}
\alias{print,MincSliceIO-method}
\alias{print,MincVolumeIO-method}
\alias{print,MincVoxelIO-method}
\alias{show-methods}
\alias{show,MincInfo-method}
\alias{show,MincSlice-method}
\alias{show,MincSliceIO-method}
\alias{show,MincVolumeIO-method}
\alias{show,MincVoxelIO-method}
\keyword{methods}

\title{Methods to Print Minc Summary Information}

\description{
The \code{print} generic function is a polymorphous function that has been overloaded
to produce summary information for the full range of \pkg{mincIO} class objects.}


\section{Methods}{
\describe{

\item{x = "MincInfo"}{Print summary information for a MincInfo object.}
\item{x = "MincSlice"}{Print summary information for a MincSlice object.}
\item{x = "MincSliceIO"}{Print summary information for a MincSliceIO object.}
\item{x = "MincVolumeIO"}{Print summary information for a MincVolumeIO object.}
\item{x = "MincVoxelIO"}{Print summary information for a MincVoxelIO object.}
}}

\usage{
\S4method{print}{MincInfo}(x, ...)

\S4method{print}{MincSlice}(x, ...)

\S4method{print}{MincSliceIO}(x, ...)

\S4method{print}{MincVolumeIO}(x, ...)

\S4method{print}{MincVoxelIO}(x, ...)
}

\arguments{
  \item{x}{The object to be summarized.}
  \item{\dots}{Other yet to be defined parameters.}
}

\value{This function is run for its side-effect.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{printMincInfo}} \code{\link{plot-methods}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# print all sorts of stuff
v305 <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")

volInfo <- readMincInfo(v305)                   #
print(volInfo)                                  # print volume info

vol <- readVolume(v305)                         #
print(vol)                                      # print volume info

# get a slice and print info
sx <- getSliceX(vol, 60)                        #
print(sx)                                       # print slice info

# read slice 80 across all frames
v4d <- file.path(volDir, "functional_4D.mnc")
sliceMatrix <- readBySlice(v4d, 80)             # load the slice matrix
print(sliceMatrix)                              # print summary info
}


