\name{civet.computeTissueVolumes}
\Rdversion{1.1}
\alias{civet.computeStxTissueVolumes}
\alias{civet.computeNativeTissueVolumes}

\title{
    Compute GM, WM, and CSF Tissue Volumes
}
\description{
    Returns a named vector of tissue volumes.
}
\usage{
civet.computeStxTissueVolumes(scanID, baseDir, civetVersion="1.1.9")

civet.computeNativeTissueVolumes(scanID, baseDir, civetVersion="1.1.9")
}

\arguments{
	\item{scanID}{A string specifying the unique scan-id (and thus sub-directory) within
			the Civet root output directory.}
	\item{baseDir}{A string specifying the Civet root output directory.  This directory 
			will, in turn, contain all of the scanIDs.}
	\item{civetVersion}{An optional string specifying the version of Civet used to create
			the output.  This is significant since filenames and directory structures may
			change across difference versions of Civet.}
}

\details{
    Actually, this function really returns the number of voxels of each tissue type contained within the final discrete Civet-produced classification volume. Now, given that Civet volumes are currently sampled using 1-mm isotropic voxels, the voxel count value should also reflect the volume in cubic millimeters. If this ever changes, we're going to have to make a minor change in this function. Please let me know if this ever happens. The native volume measurements are created by taking the stereotactic volumes and dividing each of them by the xfm-derived rescaling factor.
}

\value{A named vector containing a value for each of the 3 tissue types.}


\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}

\examples{
\dontrun{
library(RMINC)

# set Civet root path and scan-identifier
basePath <- "~/tmp/ADNI/civet/pipeOut"
scanID = "0221-M-AD"

# print gray matter volume in stereotactic space
stx_cls_vec <- civet.computeStxTissueVolumes(scanID, baseDir)
print(stx_cls_vec["gm"])

# print csf volume in native space
native_cls_vec <- civet.computeNativeTissueVolumes(scanID, baseDir)
print(native_cls_vec["csf"])
}}




