\name{minc.write.volume}
\alias{minc.write.volume}
\title{Write a buffer to a MINC2 file}
\description{Writes a buffer to a MINC2 file}
\usage{
minc.write.volume("output.mnc", "like-this-file.mnc", buffer)
}
\arguments{
  \item{output.filename}{A string containing the name of the output
    file.}
  \item{like.filename}{A string containing a filename of an existing
    MINC2 volume.}
  \item{buffer}{A 1D array containing the buffer that is to be written}
}
\details{This function writes a buffer as a MINC2 filename. The
  dimension order, step sizes, etc. will be taken from the second
  argument. Note that the data buffer, provided in the third argument,
  thus has to be of size dim1*dim2*dim3 of the like.filename for this
  function not to die a horrible death.}
\seealso{minc.get.volume, minc.apply}
