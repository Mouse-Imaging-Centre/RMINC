\name{rminc.readLinearXfmFile}
\Rdversion{1.1}
\alias{rminc.readLinearXfmFile}

\title{
    Read a Linear XFM File
}
\description{
    Read an *.xfm file, returning the contents in a data.frame.
}

\usage{
rminc.readLinearXfmFile(xfmFilename)
}

\arguments{
	\item{xfmFilename}{A string specifying the fully-qualified path to the XFM file of interest.}
}


\details{This function decodes the xfm file by making an external call to the ``xfm2param'' program within the conglomerate package.  As such, \bold{conglomerate} needs to be installed and available on the user's PATH for this function to work.
}


\value{A data.frame containing the XFM file details.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}


\examples{
\dontrun{
library(RMINC)

# set Civet root path and scan-identifier
scanID <- "0548-F-NC"
baseDir <- "~/tmp/ADNI/civet/pipeOut"

# read in the linear xfm file from Civet
xfmname <- civet.getFilenameLinearTransform(scanID, baseDir)
xfm.df <- rminc.readLinearXfmFile(xfmname)

# extract the x-scaling factor
xScale <- xfm.df["scale", "x"]

# create a global scaling factor
globalRescale <- prod(xfm.df["scale",])

}}







