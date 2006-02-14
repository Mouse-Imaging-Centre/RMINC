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
    containing group assignemnts, the correlate, or model matrix (see
    details).}
  \item{method}{A string of either "t-test", "wilcoxon", "correlation",
    or "lm" determining the method to be used.}
  \item{mask}{An optional filename containing a mask volume.}
}
\details{This function performs either a t-test, its non-parametric
  equivalent, the Wilcoxon Rank Sum test (a.k.a. the Mann-Whitney U
  test), a correlation or a linear model at every voxel. The code is entirely
  implemented in C so is much faster than calling minc.apply with the
  corresponding R functions.

  The t-test and the wilcoxon rank-sum test need subject group
  assignments as the second argument. These should be values that
  evaluate to either 0 or 1 when cast a double.

  A correlation needs a single vector of values to correlate against the
  voxel value as its second argument.

  The linear model takes the regression matrix. This is best provided by
  calling the model.matrix function here, as illustrated in the examples.
}
\value{
  \item{w}{If the method is "wilcoxon", returns a vector of length
    dim1*dim2*dim3 of W statistics.}
  \item{t}{If the method is "t-test", returns a vector of length
    dim1*dim2*dim3 of t statistics.}
  \item{r}{If the method is "correlation", returns a vector of length
    dim1*dim2*dim3 of correlation coefficients.}
  \item{lm}{If the method is "lm", returns a matrix with dim1*dim2*dim3
    rows and as many columns as there were elements in the model matrix
    of t-statistics.}
}
\seealso{minc.apply}
\examples{
\dontrun{
# get a file that could be used by glim_image
gf <- as.data.frame(read.table("filename.glim"))

# get a t-test at every voxel
t.stats <- minc.model(gf$V1, gf$V2, "t-test")

# do the non-parametric equivalent.
w.stats <- minc.model(gf$V1, gf$V2, "wilcoxon")

# run a correlation
r <- minc.model(gf$V1, gf$V3, "correlation")

# a linear model.
l <- minc.model(gf$V1, model.matrix(V1 ~ V2, gf), "lm")

# and write to file
minc.write.volume("t-test.mnc", gf$V1[1], t.stats)
minc.write.volume("w-test.mnc", gf$V1[1], w.stats)
}}
