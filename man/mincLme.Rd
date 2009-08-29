\name{mincLme}
\alias{mincLme}
\title{Linear mixed effects model at Every Voxel}
\description{Perform linear mixed effect model analysi at every
 voxel of a set of MINC files}
\usage{
out <- mincLme(files, x ~ covariate, random = ~ 1 | term, mask="mask.mnc")
}
\arguments{
  \item{formula}{The linear model formula. The left-hand term consists
    of the MINC filenames over which to compute the models at every
    voxel.}
  \item{data}{The dataframe which contains the model terms.}
  \item{random}{The random effects term to be calculated.}
  \item{subset}{Subset definition.}
  \item{mask}{Either a filename or a vector of values of the same length
    as the input files. The linear model will only be computed inside
    the mask.}
}
\details{This function computes a linear mixed effects model at every voxel of a set
  of files. The function is a close cousin to lm/lme, the key difference
  being that the left-hand side of the formula specification takes a
  series of filenames for MINC files. There is also a random term required for the analysis.}
\value{
  \item{out}{mincLme returns a mincMultiDim object which contains a series
  of columns corresponding to the terms in the linear mixed effects model. The first
  column is the F-statistic of the significance of the entire volume,
  the following columns contain the marginal t-statistics for each of
  the terms in the model}.
}
\seealso{mincWriteVolume, mincFDR, mincMean, mincSd}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLme(mincLme(files, x ~ covariate, random = ~ 1 | term)
# see what's in the results
vs
# write the results to file
mincWriteVolume(vs, "output.mnc", "Genotype+")
}}

