\name{mincRayTraceStats}
\alias{mincRayTraceStats}
\title{Create an image of a statistical peak.}
\description{Takes a voxel, an anatomical file, a statistics file
and the number of groups in the study, in order to create an image of
the statistical peak.}
\usage{
mincRayTraceStats(v, anatomy.volume, stats, groupcount,
	fdr = NULL, slice.direction = "transverse", 
	outputfile = "ray_trace_crosshair.png",
	show.pos.and.neg = FALSE, display = TRUE)
}
\arguments{
  \item{v}{A mincVoxel indicating the voxel of interest.}
  \item{anatomy.volume}{The path to the file containing the anatomical data.}
  \item{stats}{Either the path to the stats file, a mincSingleDim, or a column of a mincMultiDim.}
  \item{groupcount}{The number of groups involved in the test, to determine whether these are t-statistics of f-statistics.}
  \item{fdr}{The fdr threshold of the test.}
  \item{slice.direction }{The slice direction of the output image. This can be transverse, coronal or sagittal.}
  \item{outputfile}{The name (and path) of the outputfile.}
  \item{show.pos.and.neg}{In the case of t-statistics, when this flag is set to TRUE, the image will contain both the positive as well as the negative t-statistics.}
  \item{display}{Display the created image.}
}
\details{This function will call the ray_trace program to
create an image of a statistical peak. The anatomical slice
of the brain will be overlayed with the statistical slice and
a crosshair indicates the chosen peak.}
\seealso{mincLm, mincFDR, mincMean, mincSd}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(filenames ~ Genotype, gf)
# get the voxel at world coordinates (1,0.5,-0.5)
v <- mincGetWorldVoxel(filenames, 1, 0.5, -0.5)
# create an image of this coordinate, using the third column
# of the mincLm output.
mincRayTraceStats(v,"/some/path/anatomical.mnc", vs[,3], 2)
}}
