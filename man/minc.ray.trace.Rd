\name{minc.ray.trace}
\alias{minc.ray.trace}
\title{Call ray_trace to get a bitmap of a slice}
\description{Call ray_trace to get an image of a rendered slice}
\usage{
minc.ray.trace(volume, output="slice.rgb", size=c(400,400),
               slice=list(pos=0, wv="w", axis="z"),
               threshold=NULL,
               colourmap="-spectral",
               background=NULL,
               background.threshold=NULL,
               background.colourmap="-gray",
               display=TRUE)
}

\arguments{
  \item{volume}{The filename of a volume to render.}
  \item{output}{The output filename.}
  \item{size}{A vector of two elements specifying the output size}
  \item{slice}{A list of three elements, pos being the slice number, wv
    whether the specification is in voxel or world space, and the axis.}
  \item{threshold}{A vector of two elements containing the threshold. If
    NULL, the full range of the volume will be used.}
  \item{colourmap}{The colourmap to be used by ray_trace.}
  \item{background}{An optional filename of a background volume. Used,
    for example, to render statistical results on top of background
    anatomy.}
  \item{background.threshold}{Threshold to use for the background
    volume. If NULL the whole range will be used.}
  \item{background.colourmap}{The colourmap argument to be passed to
    ray_trace for the background image.}
  \item{display}{Boolean argument which determines whether display (from
    ImageMagick) will be called on the output.}
}
\details{
  This function provides an interface to the ray_trace command written
  by David MacDonald. As such it needs both ray_trace and make_slice to
  be on the path upon startup of R, and the bicpl library has to be
  compiled with image output enabled.

  Behaviour of minc.ray.trace varies depending on whether a background
  image is specified. If background=NULL, then the specified slice is
  rendered using the supplied (or automatically determined) threshold
  argument. If there is a background image, then the slice from the
  input volume is rendered semi-transparently on top of the background.

  Note that cropping in ray_trace is on by default, so the output image
  size will not necessarily be the same as the size argument to
  minc.ray.trace. 
}
\value{
  \item{output}{The filename of the output image is returned.}
}
\seealso{minc.apply, minc.write.volume, minc.model}
\examples{
\dontrun{
# get a file that could be used by glim_image
gf <- as.data.frame(read.table("filename.glim"))

# get a t-test at every voxel
t.stats <- minc.model(gf$V1, gf$V2, "t-test")

# and write to file
minc.write.volume("t-test.mnc", gf$V1[1], t.stats)

# create a pretty picture to include in the next Nature or Science
# article
minc.ray.trace("t-test.mnc", output="pretty.rgb", threshold=c(2.5,6),
               background="anatomy_image.mnc")
}}
