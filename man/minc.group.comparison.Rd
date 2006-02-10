\name{minc.group.comparison}
\alias{minc.group.comparison}
\title{Group Comparison at Every Voxel}
\description{Perform a group comparison between two groups at every
  voxel of a set of files}
\usage{
out <- minc.group.comparison(filenames, groupings, method=method,
mask=mask)
}
\arguments{
  \item{filenames}{A list of filenames of the MINC2 volumes across which
    to perform the group comparison.}
  \item{groupings}{A list of the same length as the number of filenames
    containing group assignemnts for each volume.}
  \item{method}{A string of either "t-test" or "wilcoxon" determining
    the method to be used.}
  \item{mask}{An optional filename containing a mask volume.}
}
\details{This function performs either a t-test or its non-parametric
  equivalent, the Wilcoxon Rank Sum test (a.k.a. the Mann-Whitney U
  test) at every voxel. The code is entirely implemented in C so is much
  faster than calling minc.apply with the corresponding R functions. The
  groupings can be factors or numbers, both of which have to evaluate to
  0 or 1. Only two groups are allowed.}
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
