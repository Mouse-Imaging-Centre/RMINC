\name{mincWriteVolume}
\alias{mincWriteVolume}
\alias{mincWriteVolume.mincSingleDim,mincWriteVolume.mincMultiDim}
\title{Write a MINC volume to file}
\usage{
  mincWriteVolume(buffer, ...)
}
\arguments{
  \item{buffer}{The data to be written to file. Usually the
  result of \link{mincLm} or some such command}

  \item{output.filename}{The filename to which to write the
  data to}

  \item{column}{Optional name of the column of a
  multidimensional MINC object to write out. By default the
  first column is used}

  \item{like.filename}{An existing MINC filename which has
  the same dimensions as the data to be written out.
  Normally this information is stored inside MINC data
  objects}

  \item{clobber}{Overwrite existing output file when set to
  TRUE, will not overwrite when set to FALSE and will
  prompt when NULL}
}
\value{
  mincLm Returns a vector of mincSingleDim class
}
\description{
  Writes a MINC volume to file
}
\details{
  This function takes numeric data, usually the results
  computed from one of the other mincFunctions, and writes
  it to file so that it can be viewed or manipulated with
  the standard MINC tools
}
\examples{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(filenames ~ Genotype, gf)
# see what's in the results
vs
# write the results to file
mincWriteVolume(vs, "output.mnc", "Genotype+")
}
\seealso{
  mincWriteVolume,mincLm,mincFDR,mincMean,mincSd
}

