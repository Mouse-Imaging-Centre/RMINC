\name{civet.functions}
\Rdversion{1.1}
\alias{civet.getFilenameClassify}
\alias{civet.getFilenameGrayMatterPve}
\alias{civet.getFilenameWhiteMatterPve}
\alias{civet.getFilenameCsfPve}
\alias{civet.getFilenameStxT1}
\alias{civet.getFilenameCerebrumMask}
\alias{civet.getFilenameSkullMask}
\alias{civet.getFilenameGrayMatterSurfaces}
\alias{civet.getFilenameWhiteMatterSurfaces}
\alias{civet.getFilenameMidSurfaces}
\alias{civet.getFilenameCorticalThickness}
\alias{civet.getFilenameMeanSurfaceCurvature}
\alias{civet.getFilenameLinearTransform}
\alias{civet.getFilenameNonlinearTransform}

\title{
    Get Selected Civet Filenames
}

\description{
    Returns either one or more Civet filenames, depending on file type.
}

\usage{
civet.getFilenameClassify(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameGrayMatterPve(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameWhiteMatterPve(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameCsfPve(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameStxT1(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameCerebrumMask(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameSkullMask(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameGrayMatterSurfaces(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameWhiteMatterSurfaces(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameMidSurfaces(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameCorticalThickness(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameMeanSurfaceCurvature(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameLinearTransform(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
civet.getFilenameNonlinearTransform(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE)
}

\arguments{
	\item{scanID}{A string specifying the unique scan-id (and thus sub-directory) within
			the Civet root output directory.}
	\item{baseDir}{A string specifying the Civet root output directory.  This directory 
			will, in turn, contain all of the scanIDs.}
	\item{civetVersion}{An optional string specifying the version of Civet used to create
			the output.  This is significant since filenames and directory structures may
			change across difference versions of Civet.}
	\item{fullPath}{A boolean specifying whether the function is to return either a 
			fully-qualified path (TRUE) or just the filename without path (FALSE).}
}

\details{
    The purpose of this function is to facilitate writing code requiring manipulation of Civet
products.  To this purpose, we have written a number of convenience functions which, given the type
of file desired and a path to the Civet output directory, are able to determine and return the 
actual filename(s).  With regard to design, the decision to include the type of file requested
in the function name rather than as a passed argument, was driven by the fact that inclusion
in the function name permits the user to easily view all file types available via an IDEs
Tab-Completion mechanism.  Thus, for example, using the OS X built-in GUI (R.app), I can type
``civet.<tab>'', and see all possible completions for all supported file types.
}

\value{Either a string or a list is returned, depending on the number of filenames returned.
		Specifically, a single filename is returned as a string, whereas multiple filenames
		are returned as named lists.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\examples{
\dontrun{
library(RMINC)

# set Civet root path and scan-identifier
basePath <- "~/tmp/ADNI/civet/pipeOut"
scanID = "0221-M-AD"

# get the name of the aggregate tissue classification volume
# ... and then read it
classifyVolname <- civet.getFilenameClassify(scanID, basePath)
classifyVol <- mincIO.readVolume(classifyVolname)

# get the left and right gray matter surface filenames and then 
# ... print the names
gmSurfName <- civet.getFilenameGrayMatterSurfaces(scanID, basePath)
print(gmSurfName$left)
print(gmSurfName$right)

# get the various transformation file filenames
lin.xfmName <- civet.getFilenameLinearTransform(scanID, basePath)
print(lin.xfmName)
nlin.xfmNames <- civet.getFilenameNonlinearTransform(scanID, basePath)
print(nlin.xfmNames$xfm)		# name of the nlin xfm file
print(nlin.xfmNames$grid)		# name of the nlin grid file
}}




