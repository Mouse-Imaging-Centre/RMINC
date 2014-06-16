\name{anatSummaries}
\alias{anatMean}
\alias{anatVar}
\alias{anatSum}
\alias{anatSd}
\title{Create descriptive statistics across a series of anat variables}

\description{This function is used to compute the mean, standard deviation,
  sum, or variance of every region specified by an anat variables.}


\usage{
anatMean(anat)

anatSd(anat)

anatSum(anat)

anatVar(anat)
}

\arguments{
  \item{filenames}{an array of atlas labels vs subject data across which to create
    the descriptive statistic.}
}

\value{
  \item{out}{The output will be a single vector containing as many
    elements as there are in the input files.}
}

\seealso{anatLm}
\examples{
\dontrun{
# read the text file describing the dataset
gf <- read.csv("control-file.csv")
# compute the mean at every voxel of all files.
means <- anatMean(gf$filenames)
}}
