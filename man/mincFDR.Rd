\name{mincFDR}
\alias{mincFDR}
\title{Compute the False Discovery Rate for a mincLm object}
\description{Takes the output of a mincLm run and computes the False
  Discovery Rate on the results.}
\usage{
qvals <- mincFDR(mincLm.object, mask="mask.mnc")
}
\arguments{
  \item{mincLm.object}{The results of a mincLm run.}
  \item{columns}{A vector of column names. By default the threshold will
    be computed for all columns; with this argument the computation can
    be limited to a subset.}
  \item{mask}{Either a filename or a numeric vector representing a mask
    - only values inside the mask will be used to compute the
    threshold.}
}
\details{This function uses the \link{qvalue} package to compute the
  False Discovery Rate threshold for the results of a \link{mincLm}
  computation. The False Discovery Rate represents the percentage of
  results expected to be a false positive.}
\value{
  \item{qvals}{mincFDR returns an object with the same number of columns
    as the input (or the subset specified by the columns argument to
    mincFDR). Each column now contains the qvalues for each voxel. Areas
    outside the mask (if a mask was specified) will be represented by a
    value of 1. The result also has an attribute called "thresholds"
    which contains the 1, 5, 10, 15, and 20 percent false discovery rate
    thresholds.}
}
\seealso{mincWriteVolume,mincLm}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(filenames ~ Genotype, gf)
# compute the False Discovery Rate
qvals <- mincFDR(vs, mask="mask.mnc")
# write the Gentoype column of the qvals to file
mincWriteVolume(qvals, "FDR-results.mnc", "Genotype+")
}}
