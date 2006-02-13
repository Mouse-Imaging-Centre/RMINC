\name{minc.model}
\alias{minc.model}
\title{Statistical model at Every Voxel}
\description{Perform one of a series of models at every voxel of a
  set of files}
\usage{
out <- minc.model(filenames, y, method=method,
mask=mask)
}
\arguments{
  \item{filenames}{A list of filenames of the MINC2 volumes across which
    to perform the statistical modelling.}
  \item{y}{A list of the same length as the number of filenames
    containing group assignemnts or the correlate.}
  \item{method}{A string of either "t-test", "wilcoxon", or
    "correlation", determining the method to be used.}
  \item{mask}{An optional filename containing a mask volume.}
}
\details{This function performs either a t-test, its non-parametric
  equivalent, the Wilcoxon Rank Sum test (a.k.a. the Mann-Whitney U
  test), or a correlation at every voxel. The code is entirely
  implemented in C so is much faster than calling minc.apply with the
  corresponding R functions.}
\seealso{minc.apply}
\examples{
\dontrun{
# get a file that could be used by glim_image
gf <- as.data.frame(read.table("filename.glim"))

# get a t-test at every voxel
t.stats <- minc.group.comparison(gf$V1, gf$V2, "t-test")

# do the non-parametric equivalent.
w.stats <- minc.group.comparison(gf$V1, gf$V2, "wilcoxon")

# and write to file
minc.write.volume("t-test.mnc", gf$V1[1], t.stats)
minc.write.volume("w-test.mnc", gf$V1[1], w.stats)
}}
