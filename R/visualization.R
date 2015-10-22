brainPlot <- function(anatomy, statistics, slice, limits=c(2,4)) {
  d <- dim(anatomy)
  as <- t(anatomy[,slice,d[3]:1])
  Sm <- melt(t(statistics[,slice,]))
  ggimage(as, fullpage=F) + 
    geom_raster(aes(x=Var2, y=Var1, fill=value), data=Sm) +
    stat_contour(aes(x=Var2, y=Var1, z=value), data=Sm) + 
    scale_fill_gradient(limits=limits, low="blue", high="red", na.value="transparent")
    #scale_fill_gradientn(colours = c(gray.colors(10), terrain.colors(10)), na.value="transparent")
}

#' A utility function to give a MINC object spatial dimensions
#'
#' Currently the plotting functions (mincPlotAnatAndStatsSlice) need an object 
#' with the correct spatial dimensions assigned. While this might in the future happen
#' at the time those objects are created, for the moment this utility function
#' works with older style RMINC objects and extract what you need.
#'
#' @param volume The input volume (from mincLm, mincGetVolume, etc.)
#' @param dimIndex The index into a multidimensional object
#'
#' @note R uses Fortran indexing, so dimension assignment is c(dim[3], dim[2], dim[1])
#' once dimensions are obtained from any libminc functions (which use C indexing)
#'
#' @return A matrix with 3 dimensions
#' @export
#'
#' @examples
#' \dontrun{
#' vol <- mincGetVolume("somefile.mnc")
#' volWithDims <- mincArray(vol)
#' 
#' vs <- mincLm(jacobians ~ genotype, gf)
#' tvol <- mincArray(vs, 6)
#' }
mincArray <- function(volume, dimIndex=1) {
  # 1d file with no dimensions (such as the output of mincGetVolume)
  if (is.null(dim(volume))) {
    outvol <- volume
  }
  # if it's already a 3d volume, just hand it back
  else if (length(dim(volume)) == 3) {
    return(volume)
  }
  # if it has two dimensions, extract just the one asked for
  else if (length(dim(volume)) == 2) {
    outvol <- volume[,dimIndex]
  }
  else {
    outvol <- volume
  }
  # we now have a 1d vector - get it's dimensions from the sizes attribute, or read it from the likeVolume
  if (! is.null(attr(volume, "sizes"))) {
    sizes <- attr(volume, "sizes")
    dim(outvol) <- c(sizes[3], sizes[2], sizes[1]) # C to Fortran dim ordering
  }
  else if (! is.null(attr(volume, "likeVolume"))) {
    sizes <- minc.dimensions.sizes(attr(volume, "likeVolume"))
    dim(outvol) <- c(sizes[3], sizes[2], sizes[1])
  }
  else {
    error("Input needs to either have dimensions set, or have a sizes or likeVolume attribute")
  }
  return(outvol)
}

getSlice <- function(volume, slice, dimension) {
  d <- dim(volume)
  if (dimension == 1) {
    outs <- volume[slice,,]
    outa <- d[3]/d[2]
  }
  else if (dimension == 2) {
    outs <- volume[,slice,]
    outa <- d[3]/d[1]
  }
  else if (dimension == 3) {
    outs <- volume[,,slice]
    outa <- d[2]/d[1]
  }
  else if (dimension > 3) {
    error("Can only handle three dimensions at the moment")
  }
  return(list(slice=outs, asp=outa))
}

mincPlotSliceSeries <- function(anatomy, statistics, dimension=2,
                                mfrow=c(4,5),
                                low=NULL, high=NULL, anatLow=NULL,
                                anatHigh=NULL, col=NULL,
                                begin=NULL, end=NULL, symmetric=F,
                                legend=NULL, plottitle=NULL) {
  opar <- par()
  nslices <- prod(mfrow)
  par(bg = "black")
  if (! is.null(legend)) { # use the layout function
    # create a layout matrix - sets up a matrix that is the same as par(mfrow=) would be,
    # and then adds on an extra column for the slice locator and the legend
    layoutMatrix <- cbind(matrix(1:nslices, nrow=mfrow[1], ncol=mfrow[2], byrow=T),
                          c(nslices+1, rep(nslices+2, mfrow[1]-1)))
    layout(layoutMatrix)
  }
  else {
    par(mfrow=mfrow)
  }
  
  # get rid of margins
  par(mar=c(0,0,0,0))
  # but keep some out margin for a title if desired
  if (!is.null(plottitle)) { 
    par(oma=c(0,0,2,0))
  }
  
  # figure out location of slices
  d <- dim(anatomy)
  if (is.null(begin)) { begin <- 1 }
  if (is.null(end)) { end <- d[dimension]-1 }
  else if (end < 0) { end <- d[dimension] + end }
  slices <- ceiling(seq(begin, end, length=nslices))
  
  # plot the actual slices
  anatRange <- getRangeFromHistogram(anatomy, anatLow, anatHigh)
  for (i in 1:nslices) {
    mincPlotAnatAndStatsSlice(anatomy, statistics, dimension, slice=slices[i],
                              low=low, high=high, anatLow=anatRange[1], 
                              anatHigh=anatRange[2], col=col, legend=NULL, 
                              symmetric=symmetric)
  }
  
  # add the slice locator and the legend if so desired
  if (!is.null(legend)) {
    if (dimension %in% c(2,3)) {
      locatorDimension = 1
      locatorSlice = ceiling(d[1]/2)
    }
    else {
      locatorDimension = 2
      locatorSlice = ceiling(d[2]/2)
    }
    mincPlotAnatAndStatsSlice(anatomy, statistics, locatorDimension, slice=locatorSlice,
                              low=low, high=high, anatLow=anatRange[1], 
                              anatHigh=anatRange[2], col=col, legend=NULL, 
                              symmetric=symmetric)
    if (dimension %in% c(1,2)) {
      abline(v=slices/d[dimension])
    }
    else {
      abline(h=slices/d[dimension])
    }
    
    # and add the colourbar 
    plot.new()
    if (symmetric==TRUE) {
      col <- colorRampPalette(c("red", "yellow"))(255)
      rcol <- colorRampPalette(c("blue", "turquoise1"))(255)
      color.legend(0.3, 0.05, 0.5, 0.45, c(high*-1, low*-1), rev(rcol), gradient="y", align="rb", col="white")
      color.legend(0.3, 0.55, 0.5, 0.95, c(low, high), col, gradient="y", align="rb", col="white")
    }
    else {
      color.legend(0.3, 0.25, 0.5, 0.75, c(low, high), col, gradient="y", align="rb", col="white")
    }
    opar <- par()
    par(xpd=T)
    text(0.85, 0.5, labels=legend, srt=90, col="white", cex=2)
  }
  if (!is.null(plottitle)) {
    mtext(plottitle, outer=T, side=3, line=-1, col="white", cex=2)
  }
  par(opar)
}

getRangeFromHistogram <- function (volume, low, high) {
  if (is.null(low)) { low <- hist(volume, plot=F)$mid[5] }
  if (is.null(high)) { high <- rev(hist(volume, plot=F)$mid)[5]}
  return(c(low, high))
}

mincTriplanarSlicePlot <- function(anatomy, statistics, slice=NULL, 
                                   layoutMatrix=NULL, ...) {
  opar <- par()
  #layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
  if (is.null(layoutMatrix)) {
    layout(matrix(c(1,2,2, 1,3,3), 2, 3, byrow=T))
  }
  par(mar=c(0,0,0,0))
  par(bg = gray.colors(255)[1])
  if (is.null(slice)) {
    slice <- ceiling(dim(anatomy)/2)
  }
  else if (length(slice) == 1) {
    slice <- rep(slice, 3)
  }
  else if (length(slice) == 3){
    # do nothing; user specified correct slices
  }
  else {
    error("not sure what to do with slice specification")
  }
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice=slice[3], dim=3, ...)
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice=slice[2], dim=2, ...)
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice=slice[1], dim=1, ...)
  par(opar)

}

mincPlotAnatAndStatsSlice <- function(anatomy, statistics, slice=NULL,
                          dimension=2, low=NULL, high=NULL,
                          anatLow=NULL, anatHigh=NULL, symmetric=F,
                          col=NULL, legend=NULL) {
  if (is.null(slice)) {
    halfdims <- ceiling(dim(anatomy)/2)
    slice <- halfdims[dimension]
  }
  if (is.null(col)) {
    if (symmetric==TRUE) {
      col <- colorRampPalette(c("red", "yellow"))(255)
      rcol <- colorRampPalette(c("blue", "turquoise1"))(255)
    }
    else {
      col <- cm.colors(255)
    }
  }
  
  anatCols = gray.colors(255, start=0.0)
    
  mincImage(anatomy, dimension, slice, axes=F, col=anatCols,
            low=anatLow, high=anatHigh)
  mincImage(statistics, dimension, slice, axes=F, add=T, col=col, underTransparent = T,
            low = low, high = high)
  if (symmetric) {
    mincImage(statistics, dimension, slice, axes=F, add=T, col=rcol, 
              underTransparent = T, reverse = T, low = low, high = high)  
  }

  if (!is.null(legend)){
    if (symmetric==TRUE) {
      color.legend(1.02, 0.05, 1.07, 0.45, c(high*-1, low*-1), rev(rcol), gradient="y", align="rb")
      color.legend(1.02, 0.55, 1.07, 0.95, c(low, high), col, gradient="y", align="rb")
    }
    else {
      color.legend(1.02, 0.25, 1.07, 0.75, c(low, high), col, gradient="y", align="rb")
    }
    opar <- par()
    par(xpd=T)
    text(1.15, 0.5, labels=legend, srt=90)
  }
}

#' Plot a slice from a MINC volume
#' 
#' Calls the \code{\link{image}} plotting function from the base R graphics with
#' some additional data munging to make it easy to make with MINC volumes
#' 
#' @param volume a MINC volume as returned by \code{\link{mincArray}}
#' @param dimension the dimension (1-3) to obtain the slice from
#' @param slice the slice number
#' @param low the low end of the range to plot. If not specified then it will be
#'   estimated based on the volume's histogram.
#' @param high the high end of the range to plot. If not specified then it will 
#'   be estimated based on the volume's histogram.
#' @param reverse whether to look only at negative numbers.
#' @param underTransparent whether to make anything below the low end of the
#'   range transparent.
#' @param ... other parameters to pass on to the \code{\link{image}} function.
#'   
#' @export
#' @examples
#' \dontrun{
#' mincImage(mincArray(anatVol), slice=100, col=gray.colors(255))
#' mincImage(mincArray(vs, 6), slice=100, col=rainbow(255), 
#'           underTransparent = T, low=2, high=6, add=T)
#' }
mincImage <- function(volume, dimension=2, slice=NULL, low=NULL, high=NULL, 
                      reverse=FALSE, underTransparent=FALSE, ...) {
  s <- getSlice(volume, slice, dimension)
  # reverse means multiply scaling by -1
  if (reverse) { m <- -1} else { m <- 1 }
  
  # determine range from histogram
  imRange <- getRangeFromHistogram(volume, low, high)
  image(scaleSlice(s$slice, imRange[1]*m, imRange[2]*m, underTransparent=underTransparent), 
        useRaster=T, asp=s$asp, ...)
}

#' Draw contour lines from a MINC volume
#'
#' @param volume the output of mincArray
#' @param dimension the dimension (from 1 to 3)
#' @param slice the slice number
#' @param ... other parameters to pass on to \code{\link{contour}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mincImage(mincArray(anatVol), slice=100, col=gray.colors(255))
#' mincContour(mincArray(anatVol), slice=100, add=T, col=rainbow(2), levels=c(1000, 1400))
#' }
mincContour <- function(volume, dimension=2, slice=NULL, ...) {
  s <- getSlice(volume, slice, dimension)
  contour(s$slice, asp=s$asp, ...)
}

scaleSlice <- function(slice, low=NULL, high=NULL, underTransparent=TRUE) {
  if (is.null(low)) {
    low <- quantile(slice, 0.5)
  }
  if (is.null(high)) {
    high <- quantile(slice, 0.7)
  }
  
  # invert if it's a negative scale
  if (high < low) {
    slice <- slice*-1
    high <- high*-1
    low <- low*-1
  }
  
  slice <- slice-low
  
  slice[slice>high] <- high
  if (underTransparent) {
    under <- NA
  }
  else { 
    under <- 0
  }
  slice[slice<0] <- under
  return(slice)
}

simpleBrainPlot <- function(anatomy, statistics, slice, dimension=2, 
                            low=NULL, high=NULL) {
  #manipulatorSetState("slice", slice)
  d <- dim(anatomy)
  aslice <- anatomy[,slice,]
  qslice <- quantile(aslice, c(0.3, 0.95))
  aslice <- aslice-qslice[1]
  aslice[aslice < 0] <- NA
  aslice[aslice>qslice[2]] <- qslice[2]
  
  sslice <- statistics[,slice,]
  sslice <- sslice-low
  sslice[sslice<=0] <- NA
  sslice[sslice>=high] <- high
  
  image(aslice, useRaster=T, col=gray.colors(100), asp=d[3]/d[1])
  image(sslice,  useRaster=T, col=rev(rainbow(100, alpha=0.7)), 
        asp=d[3]/d[1], add=T)
  #filled.contour(sslice, add=T)
}

brainLocator <- function() {
  l <- locator(n=1)
  abline(h=l$y)
  abline(v=l$x)
  l$slice <- manipulatorGetState("slice")
  return(l)
}