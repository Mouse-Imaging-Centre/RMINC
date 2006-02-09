\name{minc.apply}
\alias{minc.apply}
\title{Apply an arbitrary function to all voxels of a list of MINC
  files}
\description{Similar to apply, this allows one to apply an arbitrary R
  function to a list of MINC2 volumes.}
\usage{
average <- minc.apply(filenames, quote(mean(x)), mask=mask.filename)
}
\arguments{
  \item{filenames}{A list of filenames of MINC2 volumes.}
  \item{function.string}{A string containing the function to be
    evaluated at every voxels. See details.}
  \item{mask}{An optional filename of a mask volume.}
}
\details{This function allows for the evaluation of an arbitrary R
  function at every voxel of the input files. The function needs to take
  an argument named "x", which is an array of voxel values, and return
  an array of the same length. The examples ought to illustrate this
  somewhat better.}
\examples{
\dontrun{
# get a file that could be used by glim_image
gf <- as.data.frame(read.table("filename.glim"))

# get the mean of all files
means <- minc.apply(gf$V1, quote(mean(x)))

# run a linear model. Need to create a function first.
f <- function(x) {
  summary(lm(x ~ gf$V2))$coefficients[2,3]
}

# now apply to every voxel. This'll take a while ...
t.stats <- minc.apply(gf$V1, quote(f(x)))

# and write to file
minc.write.volume("t-test.mnc", gf$V1[1], t.stats)
}}

