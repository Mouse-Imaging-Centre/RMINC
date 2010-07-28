\name{volume.combineMaskVolumes}
\Rdversion{1.1}
\alias{volume.combineMaskVolumes}

\title{
    Combine Multiple Mask Volumes into a Single Mask Volume
}

\description{
    Given a list containing more than one label volume, combine those volumes, creating an aggregate mask volume.
}

\usage{
volume.combineMaskVolumes(vol_list)

}

\arguments{
	\item{vol_list}{A list containing more than one mask volume.  Note that all volumes must reflect the same sampling.}
}


\value{A single aggregate MincVolumeIO volume is returned.
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

# ... remove background component -- combine the rest
xl <- list(components$gm, components$wm, components$csf)
new_mask_vol <- volume.combineMaskVolumes(xl)

# display the csf+gm+wm mask
plot(new_mask_vol)

}}




