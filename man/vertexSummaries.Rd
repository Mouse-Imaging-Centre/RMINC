\name{vertexSummaries}
\alias{vertexMean}
\alias{vertexVar}
\alias{vertexSum}
\alias{vertexSd}
\title{Create descriptive statistics across a series of vertex files}

\description{This function is used to compute the mean, standard deviation,
  sum, or variance of every vertex in a set of vertex files.}


\usage{
vertexMean(filenames)

vertexSd(filenames)

vertexSum(filenames)

vertexVar(filenames)
}

\arguments{
  \item{filenames}{Filenames of the vertex volumes across which to create
    the descriptive statistic.}
}

\value{
  \item{out}{The output will be a single vector containing as many
    elements as there are vertices in the input files.}
}

\seealso{vertexLm}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# compute the mean at every voxel of all files.
means <- vertexMean(gf$filenames)
}}
