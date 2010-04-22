\name{mincIO.plot-methods}
\Rdversion{1.1}
\docType{methods}
\alias{plot-methods}
\alias{plot,MincVolumeIO,ANY-method}
\alias{plot,MincVolumeIO-method}
\alias{plot,MincSlice,ANY-method}
\alias{plot,MincSlice-method}
\keyword{methods}

\title{Methods to Plot Minc Images}

\description{
    The \code{plot} generic function is a polymorphous function that has been overloaded
    to produce 2 image plots: (1) a summary plot of an entire volume, and (2) a plot of
    a specific slice.
}

\section{Methods}{
\describe{

\item{x = "MincVolumeIO"}{A MincVolumeIO object for which to create a summary plot.  Plot consists of 21 axial slices and an
                          intensity histogram.}
\item{x = "MincSlice"}{A MincSlice object and slice to be visualized (plotted).}
}}

\usage{
\S4method{plot}{MincVolumeIO}(x, y, ...)

\S4method{plot}{MincSlice}(x, y, ...)
}

\arguments{
  \item{x}{The object to be visualized.}
  \item{y}{Not used.}
  \item{\dots}{Other yet to be defined parameters.}
}

\value{This function is run for its side-effect (plotting an image).}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\link{plotSlicePretty}} \code{\link{print-methods}}
}

\examples{
\dontrun{
# read and plot the 305 volume
v305 = "norm_avg_305_mri_minc2.mnc"

vol <-  mincIO.readVolume(v305)    # default read (grayscale)
plot(vol)                                          # ... visualize the volume
sx <-  mincIO.getSliceX(vol, 20)                   # get axial slice 20
plot(sx)                                           # ... visualize the slice

vol <-  mincIO.readVolume(v305, colorMap="rainbow")    # set display colormap on read
plot(vol)                                              # show it using "rainbow"

vol <-  mincIO.readVolume(v305, colorMap="heat.colors") #
plot(vol)                                          # now use "heat" colormap 
}}




