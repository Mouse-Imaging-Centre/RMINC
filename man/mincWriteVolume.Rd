\name{mincWriteVolume}
\alias{mincWriteVolume}
\title{Write a MINC volume to file.}
\description{Writes a MINC volume to file.}
\usage{
mincWriteVolume(buffer, "filename.mnc", "optionalColumn")
}
\arguments{
  \item{buffer}{The data to be written to file. Usually the result of
    \link{mincLm} or some such command.}
  \item{filename}{The filename to which to write the data to.}
  \item{column}{Optional name of the column of a multidimensional MINC
    object to write out. By default the first column is used.}
  \item{like.filename}{An existing MINC filename which has the same
    dimensions as the data to be written out. Normally this information
    is stored inside MINC data objects.}
}
\details{This function takes numeric data, usually the results computed
  from one of the other mincFunctions, and writes it to file so that it
  can be viewed or manipulated with the standard MINC tools.}
\seealso{mincLm,mincFDR,mincMean, mincSd}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(filenames ~ Genotype, gf)
# see what's in the results
vs
# write the results to file
mincWriteVolume(vs, "output.mnc", "Genotype+")
}}
