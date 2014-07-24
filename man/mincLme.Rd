\name{mincLme}
\alias{mincLme}
\title{Linear mixed effects model at Every Voxel}

\description{Perform linear mixed effect model analysis at every
 voxel of a set of MINC files}

\usage{
mincLme(filenames, fixed.effect, random.effect, mask=NULL)
}

\arguments{
  \item{filenames}{The filenames over which to run the ME model.}
  \item{fixed.effect}{The fixed effects term to be calculated.}
  \item{random.effect}{The random effects term to be calculated.}
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
# load lme library
library(nlme)
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# make your data frame columns directly readable by R
attach(gf)
# run a linear model relating the data in all voxels to Genotype
vs <- mincLme(files, x ~ covariate, random = ~ 1 | term)
# see what's in the results
vs
# write the results to file
mincWriteVolume(vs, "output.mnc", "Genotype+")
}}

