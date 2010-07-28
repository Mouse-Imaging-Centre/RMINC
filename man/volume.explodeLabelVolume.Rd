\name{volume.explodeLabelVolume}
\Rdversion{1.1}
\alias{volume.explodeLabelVolume}

\title{
    Explode a Label Volume into its Components
}

\description{
    Given a label volume, this function splits the volume by label and then returns a list() containing a mask volume for each of the labels.
}

\usage{
volume.explodeLabelVolume(label_vol, civetLabels=TRUE)

}

\arguments{
	\item{label_vol}{A string containing the fully-qualified path to the input label volume.}
	\item{civetLabels}{A logical variable indicating whether the label volume is using the Civet convention with regards to naming tissue type (e.g., 0=background, 1=csf, etc). If TRUE, the returned list components are named using Civet tissue types (bg, csf, gm. wm), else components are simply labelled by label number e.g. (``label_0'', ``label_2'', etc.).}
}


\value{A list is returned with each list item holding a mask volume reflecting a particular label.
}


\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\examples{
\dontrun{
library(RMINC)

# set Civet root path and scan-identifier
scanID <- "0548-F-NC"
baseDir <- "~/tmp/ADNI/civet/pipeOut"

# read in the final tissue classify volume
labelVolname <- civet.getFilenameClassify(scanID, baseDir)
label_vol <- mincIO.readVolume(labelVolname)

# explode the classify vol into its 4 components
components <- volume.explodeLabelVolume(label_vol)

# display the gray matter mask volume
plot(components$gm)

}}




