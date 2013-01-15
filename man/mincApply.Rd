\name{mincApply}
\alias{mincApply}
\alias{pMincApply}
\title{Apply arbitrary R function at every voxel}
\description{Can execute any R function at every voxel for a set of MINC
  files}
\usage{
mincApply{filenames, quote(mean(x)), mask="mask.mnc"}
pMincApply{filenames, quote(mean(x)), mask="mask.mnc", cores=4}
}
\arguments{
  \item{filenames}{The MINC files over which to apply the function. Have
    to be the same sampling.}
  \item{function string}{The function which to apply. Can only take a
    single argument, which has to be 'x'. See details and examples.}
  \item{mask}{The filename of a mask - function will only be evaluated
    inside the mask.}
  \item{maskval}{An integer specifying the value inside the mask where
    to apply the function. If left blank (the default) then anything
    above 0.5 will be considered inside the mask. This argument only
    works for mincApply, not pMincApply.}
  \item{cores}{The number of processes to split the function into. Only
    works for pMincApply, not for mincApply.}
}
\details{mincApply allows one to execute any R function at every voxel
  of a set of files. There are two variants: mincApply, which works
  inside the current R session, and pMincApply, which uses the
  multicore, doMC, and foreach packages to split the execution of the
  function across multiple cores/processors on the same machine. (Future
  work will extend pMincApply to execute across a cluster).

  Unless the function to be applied takes a single argument a wrapper
  function has to be written. This is illustrated in the examples;
  briefly, the wrapper function takes a single argument, called 'x',
  which will take on the voxel values at every voxel. The function has
  to be turned into a string; the quote function can be very handy. The
  output of this function also has to be a simple vector or scalar.

  Note that interpreted R can be very slow. Mindnumbingly slow. It
  therefore pays to write optimal functions or, whenever available, use
  the optimized equivalents. In other words and to give one example, use
  mincLm rather than applying lm, and if lm has to really be applied,
  try to use lm.fit rather than plain lm.
}
\value{
  \item{out}{The output is a matrix with the same number of rows a the
  file sizes and the same number of columns as output by the function
  that was applied. Cast into one of the MINC classes to make writing it
  out with mincWriteVolume easier.}
}
\seealso{mincWriteVolume, mincMean, mincSd, mincLm, mincAnova}
\examples{
\dontrun{
# create the wrapper function - start simple
myfunc <- function(x) {
  return(mean(x))
}
# apply and write out
out <- mincApply(gf$filenames, quote(myfunc(x)), mask="mask.mnc")
mincWriteVolume(out, "means.mnc")
}}