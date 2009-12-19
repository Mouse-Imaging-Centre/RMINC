\name{mincIO-package}
\Rdversion{1.1}
\alias{mincIO-package}
\alias{mincIO}
\docType{package}
\title{
    Read and Write MINC2 Volumes Using S4 Objects
}
\description{
    The primary purpose of this package is to provide a very easy to use
    interface for reading and writing minc2 volumes.  
}
\details{
\tabular{ll}{
Package: \tab mincIO\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-09-28\cr
License: \tab BSD\cr
LazyLoad: \tab yes\cr
}
    This package provide the user with an easy to use interface by which to
    manipulate minc2 volumes.  IO can be performed at 3 levels of granularity:
    (1) volume, (2) slice, and (3) voxel.  It is expected that most commonly
    IO would be performed at the volume level, in which either an entire 3D volume
    or one 3D frame within a 4D volume is read into an internal object containing
    a 3D array.  
    
    As volume resolution and the number of volumes to be read increases,
    so does the memory footprint of those internal structures required to hold the
    data.  In order to cope with this situation, the slice level granularity functions
    were created.  As this level of IO is performed at the slice level, many more
    volumes or frames can be loaded within a unified data structure, without running
    out of memory.  The resultant object is structured an an M x N matrix, with M 
    voxels per slice, and N slices.  The \code{\link{getSliceFromSliceArray}} function can be used 
    to convert a specific slice within the slice matrix into a slice object, for 
    visualization, or for subsequently transferring to an output volume.
}

\author{
Jim Nikelski
Maintainer: Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}
}

\keyword{ package }
