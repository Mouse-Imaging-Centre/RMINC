\name{mincIO.readMincInfo}
\Rdversion{1.1}
\alias{mincIO.readMincInfo}

\title{
    Initialize a MincInfo Object from a Valid Minc Volume
}

\description{
    This function is the primary method used for the initialization of the omnipresent MincInfo
    object.
}

\usage{
mincIO.readMincInfo(filename)
}

\arguments{
  \item{filename}{A valid Minc volume name.  Minc1 volumes are auto-magically converted to Minc2.  
                  Now, that's service.}
}

\details{
     Auto-conversion from Minc1 to Minc2.  Specifying a non-mincy volume will cause toads to rain toads 
     from the heavens.}

\value{An instantiated MincInfo object.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\seealso{
    \code{\linkS4class{MincInfo}}
}

\examples{
\dontshow{
# clear debug flag and get full volume pathname
# ... dunno why zzz.R debug setting is not passed
    R_DEBUG_mincIO <- 0      
    volDir <- system.file("packageData/volumes", package="RMINC")
}
# get information for the 305 volume
v305PET <- file.path(volDir, "average305_PET_t1_tal_lin.mnc")
vol <- mincIO.readMincInfo(v305PET)     # read the volume info
print(vol)                              # print summary info
}


