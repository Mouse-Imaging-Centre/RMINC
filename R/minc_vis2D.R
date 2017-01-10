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
    stop("Input needs to either have dimensions set, or have a sizes or likeVolume attribute")
  }
  return(outvol)
}

getSlice <- function(volume, slice, dimension) {
  d <- dim(volume)
  if(length(d) != 3) stop("volume must be 3 dimensional, you may be missing a call to mincArray")
  
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
    stop("Can only handle three dimensions at the moment")
  }
  return(list(slice=outs, asp=outa))
}

plotLocator <- function(dimension, anatomy, indicatorLevels, slices){
  d <- dim(anatomy)
  if(length(d) != 3) stop("anatomy must be 3 dimensional, you may be missing a call to mincArray")
  
  locatorDimension <- ifelse(dimension %in% 2:3, 1, 2)
  locatorSlice <- ceiling(d[locatorDimension] / 2)
  
  mincContour(anatomy, dimension=locatorDimension, 
              slice=locatorSlice, col="white", levels=indicatorLevels, 
              axes=F, drawlabels=F)
  
  if (dimension %in% 1:2) {
    abline(v=slices, col="yellow")
  }
  else {
    abline(h=slices, col="yellow")
  }
  
  invisible(NULL)
}

sliceSeriesLayout <-
  function(anatomy, dimension, mfrow, begin, end, plottitle, legend, locator){
    nslices <- prod(mfrow)
    par(bg = "black")
    if (locator || !is.null(legend)) { # use the layout function
      # create a layout matrix - sets up a matrix that is the same as par(mfrow=) would be,
      # and then adds on an extra column for the slice locator and the legend
      layoutMatrix <- cbind(matrix(1:nslices, nrow=mfrow[1], ncol=mfrow[2], byrow=T),
                            c(nslices+1, rep(nslices+2, mfrow[1]-1)))
      layout(layoutMatrix)
    }
    else { # go with standard mfrow setting
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
    if(length(d) != 3) stop("anatomy must be 3 dimensional")
  
    if (end < 0) { end <- d[dimension] + end }
    slices <- ceiling(seq(begin, end, length=nslices))
    
    return(slices)
  }

#' MINC Slice Series
#' 
#' Plot a series of slices through a minc volume
#' on a given dimension. Optionally superimpose statistics, and or include a locator
#' contour to show where slices are.
#' 
#' @param anatomy A minc array of the anatomy volume to plot
#' @param statistics optional statistics or label file to overlay on anatomy slices
#' @param dimension integer denoting which dimension to slice across
#' @param mfrow A 2 element vector of the form c(rows, columns) indicating
#' the number and position of slices to draw - slices are added by rows
#' @param low the minimum statistic to plot, taken from histogram if not supplied
#' and not \code{discreteStats}, otherwise the minimum statistic
#' @param high the maximum statistic to plot, taken from histogram if not supplied
#' and not \code{discreteStats}, otherwise the maximum statistic
#' @param anatLow the minimum anatomy intensity to plot
#' @param anatHigh the maximum antomy intensity to plot
#' @param col colours for statistics or for the anatomy if statistics are not passed
#' @param begin the first slice to plot, defaults to 1
#' @param end the last slice to plot, defaults to the last slice
#' @param symmetric whether the statistics are symmetric (such as for t-statistics)
#' @param plottitle the title of the plot if desired
#' @param legend an optional string to name the legend, indicating desire for a legend
#' (or not)
#' @param  locator whether or not to draw the locator, defaults to whether or not
#' you requested a legend
#' @param indicatorLevels numeric vector indicating where to draw slice lines on the 
#' locator, defaults to every slice
#' @param discreteStats Whether stats are discrete values and should should not have
#' their range taken from their histogram if unsupplied.
#' @details 
#' You can get a fuller tutorial on how to use the visualization tools by executing
#' the following command:
#' \code{file.show(system.file("doc/visualizationTutorial.html", package="RMINC"))}
#' 
#' On certain systems the slices are plotted with a reflected y-axis. To fix this
#' configure \code{options(RMINC_flip_image = TRUE)}  
#' @examples
#' \dontrun{
#' mincPlotSliceSeries(mincArray(anatVol),           # the anatomical volume
#'                     mincArray(vs, "tvalue-SexM"), # pull out one column of the stats
#'                     anatLow=700, anatHigh=1400,   # set anatomy thresholds
#'                     low=2.5, high=10,             # set stats thresholds
#'                     symmetric=T,                  # show separate upper and lower
#'                     begin=25, end=-25  ,          # remove slices from both sides  
#'                     legend="t-statistics")
#' }
#' @export
mincPlotSliceSeries <- 
  function(anatomy, statistics = NULL, dimension=2,
           mfrow = c(4,5),                   # layout
           low = NULL, high = NULL,          # stat thresholding
           anatLow = NULL, anatHigh = NULL,  # anatomy thresholding
           col =  heat.colors(255),
           begin = 1,                          # first slice
           end = (dim(anatomy)[dimension] - 1),# last slice 
           symmetric = FALSE,
           legend = NULL,
           locator = !is.null(legend),
           plottitle = NULL, 
           indicatorLevels = NULL,
           discreteStats = FALSE){
    
    if(length(dim(anatomy)) != 3) 
      stop("anatomy must be 3 dimensional, you may be missing a call to mincArray")
    
    plot_function <- 
      ifelse(is.null(statistics), 
             as.symbol("mincPlotSimpleSliceSeries"),
             as.symbol("mincPlotStatsSliceSeries"))
    
    dispatch_function_formals <-           #Get dispatch formal arg names
      names(formals(eval(plot_function)))  
    
    argument_list <-                       #Eval them in this environment
      lapply(dispatch_function_formals, 
             function(formal) eval(as.symbol(formal)))
    
    names(argument_list) <-                #Give the new formals their names back
      dispatch_function_formals 
                            
    do.call(eval(plot_function), argument_list)
  }


mincPlotSimpleSliceSeries <- 
  function(anatomy, dimension=2, mfrow = c(4,5),
           anatLow = NULL, anatHigh = NULL, 
           begin = 1, 
           end = (dim(anatomy)[dimension] - 1),
           plottitle = NULL,
           legend = NULL,
           locator = !is.null(legend),
           indicatorLevels = NULL){
    
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
  
    slices <- sliceSeriesLayout(anatomy, dimension, mfrow, 
                                begin, end, plottitle, legend, locator)
    
    anatRange <- getRangeFromHistogram(anatomy, anatLow, anatHigh)

    if(is.null(indicatorLevels)) indicatorLevels <- anatRange
    
    lapply(slices, function(current_slice) {
      mincImage(anatomy, dimension, slice=current_slice,
                low=anatRange[1], high=anatRange[2], 
                axes = FALSE, underTransparent = TRUE)
    })
    
    if(locator) plotLocator(dimension, anatomy, 
                            indicatorLevels, slices)
    
    if (!is.null(plottitle))
      mtext(plottitle, outer=T, side=3, line=-1, col="white", cex=2)
    
    return(invisible(NULL))
  }

mincPlotStatsSliceSeries <-
  function(anatomy, statistics, dimension=2,
           mfrow = c(4,5),                   # layout
           low = NULL, high = NULL,          # stat thresholding
           anatLow = NULL, anatHigh = NULL,  # anatomy thresholding
           col = heat.colors(255),
           begin = 1,                        # first slice
           end = dim(anatomy)[dimension] - 1,# last slice 
           symmetric = FALSE,
           legend = NULL,
           locator = !is.null(legend),
           plottitle = NULL, 
           indicatorLevels = c(900, 1200),
           discreteStats = FALSE) {
    
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    
    slices <- sliceSeriesLayout(anatomy, dimension, mfrow, 
                                begin, end, plottitle, legend, locator)
    
    anatRange <- getRangeFromHistogram(anatomy, anatLow, anatHigh)
    
    if(is.null(indicatorLevels)) indicatorLevels <- anatRange
    # force use of default colours if symmetric colours desired (for now)
    if (symmetric) {
      col <- NULL
    }
    
    # If stats are discrete don't auto-range them
    if(discreteStats){
      statRange <- c(max(min(statistics, na.rm = TRUE), low),
                     min(max(statistics, na.rm = TRUE), high))
    } else { # otherwise do
      statRange <- getRangeFromHistogram(statistics, low, high)
    }

    for (i in 1:length(slices)) {
      mincPlotAnatAndStatsSlice(anatomy, statistics, dimension, slice=slices[i],
                                low=statRange[1], high=statRange[2], anatLow=anatRange[1], 
                                anatHigh=anatRange[2], col=col, legend=NULL, 
                                symmetric=symmetric)
    }
    
    #Add the plot locator if desired
    if(locator) plotLocator(dimension, anatomy, 
                            indicatorLevels, slices)
    
    #Add a legend for the statistics if desired
    if(!is.null(legend)){
      plot.new()
      if (symmetric==TRUE) {
        col <- colorRampPalette(c("red", "yellow"))(255)
        rcol <- colorRampPalette(c("blue", "turquoise1"))(255)
        plotrix::color.legend(0.3, 0.05, 0.5, 0.45, c(high*-1, low*-1), 
                              rev(rcol), gradient="y", align="rb", col="white")
        plotrix::color.legend(0.3, 0.55, 0.5, 0.95, c(low, high), 
                              col, gradient="y", align="rb", col="white")
      }
      else {
        plotrix::color.legend(0.3, 0.25, 0.5, 0.75, c(low, high), 
                               col, gradient="y", align="rb", col="white")
      }
      
      text(0.85, 0.5, labels=legend, srt=90, col="white", cex=2)
    }
    
    #Add a title if desired
    if (!is.null(plottitle)) {
      mtext(plottitle, outer=T, side=3, line=-1, col="white", cex=2)
    }
    
    return(invisible(NULL))
  }

getRangeFromHistogram <- function (volume, low = NULL, high = NULL) {
  
  if(is.null(low) || is.null(high)) hist_midpoints <- hist(volume, plot=F)$mids
  if (is.null(low)) { low <- hist_midpoints[5]}
  if (is.null(high)) { high <- rev(hist_midpoints)[5]}
  
  return(c(low, high))
}

#' Plot Slice Along Each Axis
#' 
#' Show a slice from each axis of minc volume
#' 
#' @param anatomy a \link{mincArray} object containing the source anatomy
#' @param statistics a \link{mincArray} object containing a statistic to overlay
#' @param slice 3-component vector indicating which slice along each axis
#' @param layoutMatrix A matrix describing the layout for the plots typically
#' produced by \link{layout}
#' @param ... extra parameters to be passed to \link{mincPlotAnatAndStatsSlice}
#' @return invisible NULL 
#' @export 
mincTriplanarSlicePlot <- function(anatomy, statistics, slice=NULL, 
                                   layoutMatrix=NULL, ...) {
  opar <- par(no.readonly = TRUE)
  
  if(length(dim(anatomy)) != 3) 
    stop("anatomy must be 3 dimensional, you may be missing a call to mincArray")
  
  #layout(matrix(c(1,1,2,3), 2, 2, byrow=T))
  if (is.null(layoutMatrix)) {
    layout(matrix(c(1,2,2, 1,3,3), 2, 3, byrow=T))
  }
  par(mar=c(0,0,0,0))
  par(bg = gray.colors(255, start=0)[1])
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
    stop("not sure what to do with slice specification")
  }
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice = slice[3], dimension = 3, ...)
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice = slice[2], dimension = 2, ...)
  mincPlotAnatAndStatsSlice(anatomy, statistics, slice = slice[1], dimension = 1, ...)
  par(opar)

  invisible(NULL)
}

# note - works, but is extremely slow. In profiling it appears to spend almost all
# its time in "save", so I think the function has to be rewritten in such a way as to avoid
# constant saving of large data to disk. Somehow.
# mincPlotAnatAndStatsSliceManipulator <- function(anatomy, statistics,
#                                                  slice=NULL,
#                                                  dimension=2,
#                                                  low=NULL,
#                                                  high=NULL,
#                                                  anatLow=NULL,
#                                                  anatHigh=NULL,
#                                                  legend=NULL,
#                                                  symmetric= FALSE,
#                                                  discreteStats = FALSE) {
#   d <- dim(anatomy)
#   maxStats <- max(abs(statistics))
#   maxAnat <- max(anatomy)
#   manipulate::manipulate(
#     mincPlotAnatAndStatsSlice(anatomy, statistics,
#                               slice=xSlice,
#                               dimension=dimension,
#                               low=xLow,
#                               high=xHigh,
#                               anatLow=xanatLow,
#                               anatHigh=xanatHigh,
#                               legend=legend,
#                               symmetric=xSymmetric),
#     xSlice = slider(1, d[dimension], initial=slice, label="Slice"),
#     xLow =slider(0, maxStats, initial=low, label="lower(statistics)"),
#     xHigh = slider(0, maxStats, initial=high, label="upper(statistics)"),
#     xanatLow = slider(0, maxAnat, initial=anatLow, label="lower(anatomy)"),
#     xanatHigh = slider(0, maxAnat, initial=anatHigh, label="upper(anatomy)"),
#     xSymmetric = checkbox(initial=symmetric, label="symmetrix")
#   )
# }

#' Anatomy and Statistics Slice
#' 
#' Plot a given slice through a MINC volume and superimpose statistics
#' on the slice.
#' 
#' @param anatomy A minc array of the anatomy volume to plot
#' @param statistics optional statistics or label file to overlay on anatomy slices
#' @param slice the voxel index of the slice of interest
#' @param dimension integer denoting which dimension to slice across
#' @param low the minimum statistic to plot
#' @param high the maximum statistic to plot
#' @param anatLow the minimum anatomy intensity to plot
#' @param anatHigh the maximum antomy intensity to plot
#' @param symmetric whether the statistics are symmetric (such as for t-statistics)
#' @param col colours for statistics
#' @param rcol colours for negative statistics if using a symmetric statistic
#' @param legend an optional string to name the legend, indicating desire for a legend
#' (or not)
#' @return invisible NULL
#' @export
mincPlotAnatAndStatsSlice <- function(anatomy, statistics, slice=NULL,
                          dimension=2, 
                          low=min(statistics, na.rm = TRUE), 
                          high=max(statistics, na.rm = TRUE), 
                          anatLow=min(anatomy, na.rm = TRUE), 
                          anatHigh=max(anatomy, na.rm = TRUE), 
                          symmetric=FALSE,
                          col=NULL, rcol=NULL, legend=NULL) {
  
  if(length(dim(anatomy)) != 3) 
    stop("anatomy must be 3 dimensional, you may be missing a call to mincArray")
  
  if (is.null(slice)) {
    halfdims <- ceiling(dim(anatomy)/2)
    slice <- halfdims[dimension]
  }
  if (is.null(col)) {
    if (symmetric==TRUE) {
      col <- colorRampPalette(c("red", "yellow"))(255)
    } else {
      col <- rainbow(255)
    }
  }
  if (is.null(rcol) && symmetric) {
    rcol <- colorRampPalette(c("blue", "turquoise1"))(255)
  }
  
  anatCols = gray.colors(255, start=0.0)
    
  mincImage(anatomy, dimension, slice, axes = FALSE, col=anatCols,
            low=anatLow, high=anatHigh)
  
  mincImage(statistics, dimension, slice, axes = FALSE, add = TRUE, col=col, underTransparent = TRUE,
            low = low, high = high)
  
  if (symmetric) {
    mincImage(statistics, dimension, slice, axes = FALSE, add = TRUE, col=rcol, 
              underTransparent = TRUE, reverse = TRUE, low = low, high = high)  
  }

  if (!is.null(legend)){
    # get the dimensions of the current plot
    plotdims <- par("usr")
    if (symmetric==TRUE) {
      plotrix::color.legend(0.97 * plotdims[2], 
                             0.05 * plotdims[4], 
                             0.99 * plotdims[2], 
                             0.45 * plotdims[4], 
                             c(high*-1, low*-1), rev(rcol), gradient="y", align="rb")
      plotrix::color.legend(0.97 * plotdims[2], 
                             0.55 * plotdims[4], 
                             0.99 * plotdims[2], 
                             0.95 * plotdims[4], 
                             c(low, high), col, gradient="y", align="rb")
      text(1.10, 0.5, labels=legend, srt=90)
    }
    else {
      plotrix::color.legend(0.97 * plotdims[2], 
                             0.25 * plotdims[4], 
                             0.99 * plotdims[2], 
                             0.75 * plotdims[4], 
                             c(low, high), col, gradient="y", align="rb")
      text(1.05, 0.5, labels=legend, srt=90)
    }
    opar <- par(no.readonly = TRUE) #I think this gets lost anyway
    par(xpd=T)
  }
  
  invisible(NULL)
}

#' Plot a slice from a MINC volume
#' 
#' Calls the \code{\link{image}} plotting function from the base R graphics with
#' some additional data munging to make it easy to work with MINC slices
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
#' @param col a colour palette AKA look-up-table to colourize the slice intensities
#' @param add whether to add the slice to the current plot device or open a new 
#' one before plotting 
#' @param ... other parameters to pass on to the \code{\link{image}} function.
#' @details 
#' You can get a fuller tutorial on how to use the visualization tools by executing
#' the following command:
#' \code{file.show(system.file("doc/visualizationTutorial.html", package="RMINC"))}
#' 
#' On certain systems the slices are plotted with a reflected y-axis. To fix this
#' configure \code{options(RMINC_flip_image = TRUE)}   
#' @export
#' @examples
#' \dontrun{
#' mincImage(mincArray(anatVol), slice=100, col=gray.colors(255))
#' mincImage(mincArray(vs, 6), slice=100, col=rainbow(255), 
#'           underTransparent = T, low=2, high=6, add=T)
#' }
mincImage <- function(volume, dimension=2, slice=NULL, 
                      low=min(volume, na.rm = TRUE), 
                      high=max(volume, na.rm = TRUE), 
                      reverse=FALSE, underTransparent=FALSE, 
                      col = gray.colors(255),
                      add = FALSE,
                      ...) {
  if(length(dim(volume)) != 3) 
    stop("volume must be 3 dimensional, you may be missing a call to mincArray")
  
  s <- getSlice(volume, slice, dimension)
  # reverse means multiply scaling by -1
  if (reverse) { m <- -1 } else { m <- 1 }
  
  imRange <- getRangeFromHistogram(volume, low, high)
  s$slice <- scaleSlice(s$slice, imRange[1]*m, imRange[2]*m, 
                        underTransparent=underTransparent)
  
  sliceDims <- dim(s$slice)
  
  if(!add){
    plot.special <- #override plot.defaults defaults
      function(..., xlab = "", ylab = "", asp = 1, 
               xaxs = "i", yaxs = "i")
        plot.default(..., xlab = xlab, ylab = ylab, asp = asp,
                     xaxs = xaxs, yaxs = yaxs)
    
    plot.special(
      sliceDims[1], sliceDims[2], #dummy plot values
      type = "n", #don't plot them
      xlim = c(0, sliceDims[1]),
      ylim = c(0, sliceDims[2]),
      ...)
  }
  
  if(!all(is.na(s$slice))){
    colourDepth <- length(col)

    paletteScaledSlice <- scaleSliceToPalette(s$slice, low, high, col)
    
    colourizedSlice <- col[paletteScaledSlice]
    dim(colourizedSlice) <- sliceDims
    
    flip_option <- getOption("RMINC_flip_image", FALSE)
    if(flip_option)
      colourizedSlice <- colourizedSlice[,sliceDims[2]:1]
      
    colourizedSlice <- t(colourizedSlice) #transpose for raster plotting
    
    rasterImage(colourizedSlice, 
                xleft = 0, xright = sliceDims[1],
                ytop = 0, ybottom = sliceDims[2])
  }
  
  return(invisible(NULL))
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
  
  slice[slice >= high] <- high
  slice <- slice - low
  
  if (underTransparent) {
    under <- NA
  }
  else { 
    under <- 0
  }
  
  slice[slice <= 0] <- under
  return(slice)
}

scaleSliceToPalette <- 
  function(slice, low, high, palette){
    dims <- dim(slice)
    slice <- #Scale to 0-1
      slice / abs(high - low)
    
    maxima <- which(slice == 1)
    slice[maxima] <- slice[maxima] - .Machine$double.eps
    slice <- slice * length(palette)
    slice <- floor(slice) + 1
    dim(slice) <- dims
    
    return(slice)
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
  
  if(length(dim(volume)) != 3) 
    stop("volume must be 3 dimensional, you may be missing a call to mincArray")
  
  s <- getSlice(volume, slice, dimension)
  sliceDims <- dim(s$slice)
  contour(1:sliceDims[1], 1:sliceDims[2], s$slice, asp=1, ...)
}

# simpleBrainPlot <- function(anatomy, statistics, slice, dimension=2, 
#                             low=NULL, high=NULL) {
#   #manipulatorSetState("slice", slice)
#   d <- dim(anatomy)
#   aslice <- anatomy[,slice,]
#   qslice <- quantile(aslice, c(0.3, 0.95))
#   aslice <- aslice-qslice[1]
#   aslice[aslice < 0] <- NA
#   aslice[aslice>qslice[2]] <- qslice[2]
#   
#   sslice <- statistics[,slice,]
#   sslice <- sslice-low
#   sslice[sslice<=0] <- NA
#   sslice[sslice>=high] <- high
#   
#   image(aslice, useRaster=T, col=gray.colors(100), asp=d[3]/d[1])
#   image(sslice,  useRaster=T, col=rev(rainbow(100, alpha=0.7)), 
#         asp=d[3]/d[1], add=T)
#   #filled.contour(sslice, add=T)
# }
# 
# brainLocator <- function() {
#   l <- locator(n=1)
#   abline(h=l$y)
#   abline(v=l$x)
#   l$slice <- manipulate::manipulatorGetState("slice")
#   return(l)
# }
# 
# my.get.coord <- function() {
#   par(mfg = c(1,1)) #locator() shall be relative to the first plot out
#   # of the eight plots totally
#   my.loc <-locator(1) #location, not in inches
#   my.plot.region <- par("usr") #extremes of plotting region
#   #(in plot units, not inches)
#   my.plot.region.x <- my.plot.region[2] - my.plot.region[1]
#   my.plot.region.y <- my.plot.region[4] - my.plot.region[3]
#   my.loc.inch.x <- (my.loc$x + 0.5)/my.plot.region.x * (par("pin")[1]) 
#   #par("pin") #current plot dimension in inches
#   #relative to the plotting-region bottom left corner, not the axis c(0,0) point
#   my.loc.inch.y <- (my.loc$y + 0.5)/my.plot.region.y * (par("pin")[2])
#   
#   ## search the plot we are in with locator(1)
#   my.plot.inch.x <- par("pin")[1] + par("mai")[2] + par("mai")[4] #plot.x + left & right margin
#   par("fin")[1]
#   my.plot.inch.y <- par("pin")[2] + par("mai")[1] + par("mai")[3] #plot.y + bottom & top margin
#   par("fin")[2]
#   
#   pos.rel.x <- (my.loc.inch.x / par("fin")[1] - floor(my.loc.inch.x / 
#                                                         par("fin")[1])) *
#     par("fin")[1] / par("pin")[1] * (par("usr")[2] - par("usr")[1]) - 0.5
#   #inches from left bottom corner in target plot region (c(0,0)
#   # is plot-region bottom-left corner, not the axis c(0,0) point
#   pos.rel.y <- (my.loc.inch.y / par("fin")[2] - floor(my.loc.inch.y / 
#                                                         par("fin")[2])) *
#     par("fin")[2] / par("pin")[2] * (par("usr")[4] - par("usr")[3]) - 0.5
#   #inches from left bottom corner in target plot
#   
#   fig.coord.x <- ceiling(my.loc.inch.x / par("fin")[1])
#   fig.coord.y <- 1 +(-1) *ceiling(my.loc.inch.y / par("fin")[2])
#   # cat("figure-coord x: ", fig.coord.x,"\n")
#   # cat("figure-coord y: ", fig.coord.y,"\n")
#   cat("we are in figure: ", fig.coord.y * nc + fig.coord.x, "\n")
#   cat("coordinates of the identified point x: ", pos.rel.x,"\n")
#   cat("coordinates of the identified point y: ", pos.rel.y,"\n")
# }

#' Call ray_trace to get an image of a rendered slice
#' 
#' This function provides an interface to the ray_trace command written by
#' David MacDonald. As such it needs both ray_trace and make_slice to be on the
#' path upon startup of R, and the bicpl library has to be compiled with image
#' output enabled.
#' 
#' Behaviour of minc.ray.trace varies depending on whether a background image
#' is specified. If background=NULL, then the specified slice is rendered using
#' the supplied (or automatically determined) threshold argument. If there is a
#' background image, then the slice from the input volume is rendered
#' semi-transparently on top of the background.
#' 
#' Note that cropping in ray_trace is on by default, so the output image size
#' will not necessarily be the same as the size argument to minc.ray.trace.
#' 
#' @param volume The filename of a volume to render.
#' @param output The output filename.
#' @param size A vector of two elements specifying the output size
#' @param slice A list of three elements, pos being the slice number, wv
#' whether the specification is in voxel or world space, and the axis.
#' @param threshold A vector of two elements containing the threshold. If NULL,
#' the full range of the volume will be used.
#' @param colourmap The colourmap to be used by ray_trace.
#' @param background An optional filename of a background volume. Used, for
#' example, to render statistical results on top of background anatomy.
#' @param background.threshold Threshold to use for the background volume. If
#' NULL the whole range will be used.
#' @param background.colourmap The colourmap argument to be passed to ray_trace
#' for the background image.
#' @param display Boolean argument which determines whether display (from
#' ImageMagick) will be called on the output.
#' @return \item{output}{The filename of the output image is returned.}
#' @export
minc.ray.trace <- function(volume, output="slice.rgb", size=c(400,400),
                           slice=list(pos=0, wv="w", axis="z"),
                           threshold=NULL,
                           colourmap="-spectral",
                           background=NULL,
                           background.threshold=NULL,
                           background.colourmap="-gray",
                           display=TRUE) {
  # create the slice obj
  slice.obj.name <- "/tmp/slice.obj"
  system(paste("make_slice", volume, slice.obj.name, slice$axis,
               slice$wv, slice$pos, sep=" "))
  
  # get the threshold if necessary
  if (is.null(threshold)) {
    vol <- mincGetVolume(volume)
    threshold <- range(vol)
    rm(vol)
  }
  
  # get the background threshold if necessary
  if (!is.null(background) && is.null(background.threshold)) {
    vol <- mincGetVolume(background)
    background.threshold <- range(vol)
    rm(vol)
  }
  
  # call ray_trace
  position <- ""
  if (slice$axis == "y") {
    position <- "-back"
  }
  
  if (is.null(background)) {
    system(paste("ray_trace -output", output, colourmap, threshold[1],
                 threshold[2], volume, "0 1", slice.obj.name,
                 "-bg black -crop -size", size[1], size[2], position))
  } else {
    system(paste("ray_trace -output", output, background.colourmap,
                 background.threshold[1], background.threshold[2],
                 background, "0 1 -under transparent", colourmap,
                 threshold[1], threshold[2], volume, "0 0.5", slice.obj.name,
                 "-bg black -crop -size", size[1], size[2], position))
  }
  if (display) {
    system(paste("display", output))
  }
  return(output)
}

#' Create an image of a statistical peak.
#' 
#' Takes a voxel, an anatomical file, and a statistics file to create an image
#' of the statistical peak.
#' 
#' This function will call the ray_trace program to create an image of a
#' statistical peak. The anatomical slice of the brain will be overlayed with
#' the statistical slice and a crosshair indicates the chosen peak.
#' 
#' @param v A mincVoxel indicating the voxel of interest.
#' @param anatomy.volume The path to the file containing the anatomical data.
#' @param statsbuffer Either the path to the stats file, a mincSingleDim, or a
#' mincMultiDim.
#' @param column If a mincMultiDim is specified, column will indicate which
#' column of the data to use.
#' @param like.filename If a column of a mincMultiDim is explicity passed
#' through, a like file is needed to determine the dimensions of the stats
#' buffer.
#' @param mask If a mask is specified, the stats outside the mask will not be
#' displayed in output image.
#' @param image.min Specify the minimum image intensity.
#' @param image.max Specify the maximum image intensity.
#' @param output.width Specify the width of the output image.
#' @param output.height Specify the height of the output image.
#' @param place.inset Boolean indicating whether or not to place a 3D brain
#' inset.
#' @param inset Path to the object file (.obj) containing the surface of the
#' brain.
#' @param stats.largest.pos Specify the maximum stats value.
#' @param stats.largest.neg Specify the minimum stats value.
#' @param caption Specify the caption for the colourbar. If spaces occur in the
#' caption use sometime along the line caption="\"Captoin with spaces\"".
#' @param fdr Specify the statistical significance threshold.
#' @param slice.direction The slice direction of the output image. This can be
#' transverse, coronal or sagittal.
#' @param outputfile The name (and path) of the outputfile.
#' @param show.pos.and.neg In the case of t-statistics, when this flag is set
#' to TRUE, the image will contain both the positive as well as the negative
#' t-statistics.
#' @param display Display the created image.
#' @param clobber Overwrite existing output file when set to TRUE, will not
#' overwrite when set to FALSE and will prompt when NULL.
#' @param tmpdir Specify a directory for temporary files.
#' @seealso mincLm, mincFDR, mincMean, mincSd
#' @examples
#' 
#' \dontrun{
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(filenames ~ Genotype, gf)
#' # get the voxel at world coordinates (1,0.5,-0.5)
#' v <- mincGetWorldVoxel(filenames, 1, 0.5, -0.5)
#' # create an image of this coordinate, using the third column
#' # of the mincLm output.
#' mincRayTraceStats(v,"/some/path/anatomical.mnc", vs[,3], like.filename = "like-this-file.mnc")
#' # in this particular case, a like file is stored with the vs object and
#' # can be retrieved using:
#' mincRayTraceStats(v,"/some/path/anatomical.mnc", vs[,3], like.filename = attr(vs, "likeVolume"))
#' # or
#' mincRayTraceStats(v,"/some/path/anatomical.mnc", vs, column = 3)
#' }
#' @export
mincRayTraceStats <- function(v, anatomy.volume, 
                              statsbuffer, column=1, like.filename=NULL, 
                              mask=NULL, image.min=-1000, image.max=4000,
                              output.width=800, output.height=800, 
                              place.inset=FALSE, inset=NULL,
                              stats.largest.pos=NULL, stats.largest.neg=NULL,
                              caption="t-statistic",
                              fdr=NULL, slice.direction="transverse",
                              outputfile="ray_trace_crosshair.png", 
                              show.pos.and.neg=FALSE, display=TRUE,
                              clobber=NULL, tmpdir="/tmp"){
  #check whether ray_trace_crosshair is installed
  lasterr <- try(system("ray_trace_crosshair", ignore.stderr = TRUE), 
                 silent=TRUE)
  if(lasterr == 32512){
    stop("ray_trace_crosshair must be installed for mincRayTraceStats to work.")
  }
  
  if(file.exists(outputfile) && is.null(clobber)){
    answer <- readline("Warning: the outputfile already exists, continue? (y/n) ")
    if(substr(answer, 1, 1) == "n")
      stop("Output file exists, specify clobber, or change the output file name.")
  }
  else if(file.exists(outputfile) && !clobber){
    stop("Output file exists, specify clobber, or change the output file name.")
  }
  
  ### VOXEL
  #check whether the first argument is a mincVoxel:
  if(class(v)[1] != "mincVoxel"){
    stop("Please specify a mincVoxel")
  }
  
  #create a system call for ray_trace_crosshair.pl
  systemcall <- array()	
  systemcall[1] <- "ray_trace_crosshair"
  systemcall[2] <- "-x"
  systemcall[3] <- attr(v,"worldCoord")[1]
  systemcall[4] <- "-y"
  systemcall[5] <- attr(v,"worldCoord")[2]
  systemcall[6] <- "-z"
  systemcall[7] <- attr(v,"worldCoord")[3]
  systemcall[8] <- "-caption"
  systemcall[9] <- caption
  
  ### ANATOMY VOLUME
  i <- 10;
  systemcall[i] <- "-final-nlin"
  i <- i + 1
  if(class(anatomy.volume)[1] == "character"){
    systemcall[i] <- anatomy.volume
    i <- i + 1
  }
  else if(class(anatomy.volume)[1] == "mincSingleDim"){
    systemcall[i] <- attr(anatomy.volume, "filename")
    i <- i + 1
  }
  else{
    stop("Please specify the path to the anatomy volume, or a mincSingleDim containing the anatomy volume.")
  }
  
  ### STATS VOLUME
  systemcall[i] <- "-jacobian"
  i <- i + 1
  if(class(statsbuffer)[1] == "character"){
    systemcall[i] <- statsbuffer
    i <- i + 1
  }
  
  if(class(statsbuffer)[1] == "numeric"){
    if (is.na(file.info(as.character(like.filename))$size)){
      stop(c("File ", like.filename, " cannot be found.\n"))
    }
    #write buffer to file
    mincWriteVolume.default(statsbuffer, 
                            paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep=""), 
                            like.filename)
    
    systemcall[i] <- paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep="")
    i <- i + 1
  }
  
  if(class(statsbuffer)[1] == "mincMultiDim"){
    if (is.null(like.filename)){
      like.filename <- attr(statsbuffer, "likeVolume")
    }
    if (is.na(file.info(like.filename)$size)){
      stop(c("File ", like.filename, " cannot be found.\n"))
    }
    
    #write buffer to file
    mincWriteVolume.default(statsbuffer[,column], 
                            paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep=""),
                            like.filename)
    
    systemcall[i] <- paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep="")
    i <- i + 1
  }
  
  if(class(mask) != "NULL"){
    path.to.mask = "init"
    if(class(mask)[1] == "character"){
      path.to.mask = mask
    }
    else if(class(mask)[1] == "mincSingleDim"){
      path.to.mask <- attr(mask, "filename")
    }
    system(paste("mv", 
                 paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep=""),
                 paste(tmpdir, "/R-wrapper-ray-trace-stats-full.mnc", sep="")))
    system(paste("mincmask", 
                 paste(tmpdir, "/R-wrapper-ray-trace-stats-full.mnc", sep=""),
                 path.to.mask,
                 paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep="")))
    system(paste("rm -f",
                 paste(tmpdir, "/R-wrapper-ray-trace-stats-full.mnc", sep="")))
  }
  
  ### IMAGE INTENSITY EXTREMA
  systemcall[i] <- "-image-min"
  i <- i + 1
  systemcall[i] <- image.min
  i <- i + 1
  systemcall[i] <- "-image-max"
  i <- i + 1
  systemcall[i] <- image.max
  i <- i + 1
  
  ### OUTPUT IMAGE DIMENSIONS
  systemcall[i] <- "-image-width"
  i <- i + 1
  systemcall[i] <- output.width
  i <- i + 1
  systemcall[i] <- "-image-height"
  i <- i + 1
  systemcall[i] <- output.height
  i <- i + 1
  
  ### INSET
  if(place.inset == FALSE){
    systemcall[i] <- "-no-place-inset"
    i <- i + 1
  }
  else{
    systemcall[i] <- "-place-inset"
    i <- i + 1
  }
  
  if(! is.null(inset)){
    systemcall[i] <- "-brainsurface"
    i <- i + 1
    systemcall[i] <- inset
    i <- i + 1
  }
  
  ### STATS EXTREMA
  if(! is.null(stats.largest.pos)){
    systemcall[i] <- "-positive-max"
    i <- i + 1
    systemcall[i] <- stats.largest.pos
    i <- i + 1
  }
  
  if(! is.null(stats.largest.neg)){
    systemcall[i] <- "-negative-min"
    i <- i + 1
    systemcall[i] <- stats.largest.neg
    i <- i + 1
  }
  
  
  ### FDR
  if(class(fdr)[1] == "numeric"){
    systemcall[i] <- "-fdr"
    i <- i + 1
    systemcall[i] <- fdr
    i <- i + 1
  }
  
  ### OUTPUTFILE
  systemcall[i] <- "-outputfile"
  i <- i + 1
  systemcall[i] <- outputfile
  i <- i + 1
  
  ### SLICE DIRECTION
  systemcall[i] <- "-slicedirection"
  i <- i + 1
  systemcall[i] <- slice.direction
  i <- i + 1
  
  ### SHOW POSITIVE AND NEGATIVE
  if(show.pos.and.neg == FALSE){
    systemcall[i] <- "-no-show-pos-and-neg"
    i <- i + 1
  }
  else{
    systemcall[i] <- "-show-pos-and-neg"
    i <- i + 1
  }
  
  ### #### ###
  systemcall[i] <- "-remove-temp"
  i <- i + 1
  system(paste(systemcall, collapse = " " ))
  
  if(display){
    system(paste("display", outputfile, "&"))
  }
  
  #remove the file written to disk
  system(paste("rm -f",
               paste(tmpdir, "/R-wrapper-ray-trace-stats.mnc", sep="")))
  
}


#' Plotting of peaks
#'
#' Plots a slice containing a peak. Optionally plots a graph of
#' that peak alongside.
#'
#' @param peak a row from \code{\link{mincFindPeaks}}
#' @param anatomy a mincArray for the underlying anatomy
#' @param statistics a mincArray for the stats volume
#' @param dim the dimension (1:3)
#' @param crossCol the colour for the cross-hair
#' @param crossSize the size (cex) of the cross-hair
#' @param plotFunction a function which will produce a graph
#' @param ... other details to pass on to \code{\link{mincPlotAnatAndStatsSlice}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' peaks <- mincFindPeaks(-log10(qvs), "Neonatal:time.to.sac", "pos",
#' posThreshold=1.3, minDistance=1)
#' p <- function(peak) {
#'   gfTiming$voxel <- mincGetWorldVoxel(gfTiming$reljacobians02,
#'                                       peak["x"], peak["y"], peak["z"])
#'   qplot(time.to.sac, exp(voxel), data=gfTiming, colour=Neonatal,
#'         geom="boxplot") + theme_classic()
#' }
#' mincPlotPeak(peaks[1,], anatVol, -log10(mincArray(qvs, "Neonatal:time.to.sac")), 
#'              anatLow=700, anatHigh=1400, low=1, high=4, col=heat.colors(244), 
#'              crossCol = "blue", crossSize = 3, plotFunction = p)
#' }
mincPlotPeak <- function(peak, anatomy, statistics, dim=2, 
                         crossCol = "green", crossSize=4,
                         plotFunction=NULL, ...) {
  
  # get rid of margins
  oldpar <- par(mar=c(0,0,0,0))
  # if we have a plotFunction, divide plot into two
  if (!is.null(plotFunction)) {
    oldpar <- par(mfrow=c(1,2))
  }
  # at it's most basic, show the slice with a cross-hair
  mincPlotAnatAndStatsSlice(anatomy, statistics, dimension = dim, 
                            slice = as.numeric(peak[dim]), ...)
  if (dim == 1) {
    points(peak$d2, peak$d3, col=crossCol, pch="+", cex=crossSize)
  } else if (dim == 2) {
    points(peak$d1, peak$d3, col=crossCol, pch="+", cex=crossSize)
  }
  else {
    points(peak$d1, peak$d2, col=crossCol, pch="+", cex=crossSize)
  }
  # add the plot 
  if (!is.null(plotFunction)) {
    pp <- plotFunction(peak)
    # if the plot is grid (ggplot2), then need to play with viewports
    if( any(class(pp) %in% "gg") ) {
      plot.new()
      vps <- baseViewports()
      pushViewport(vps$figure)
      vp1 <- plotViewport(c(1,1,1,1))
      print(pp, vp=vp1)
      popViewport()
    }
  }
  par(oldpar)
}

#' A tool that returns a color function/palette from color lookup files
#'
#' @param lookup_table Either a path to the lookup table file, or the table itself
#' @param alpha A transparency value between 0 and 1 (inclusive)
#' 
#' @note For now, the first column of the lookup table (interval where the colour occurs) is ignored; i.e. equal spacing intervals are assumed
#' 
#' @return A function that takes in an integer specifying the number of colours required, and returns a vector of colours interpolated from the input colour lookup file
#' @export
#'
#' @examples
#' \dontrun{
#' spectral.colors <- mincColorLookupToColorRampPalette("/micehome/jlerch/luts/spectral")
#' spectral.colors(100)
#' }
mincColorLookupToColorRampPalette <- function(lookup_table="/micehome/jlerch/luts/spectral", alpha=1) {
  tryCatch({
    lut <- read.table(file=lookup_table, header=FALSE, sep = " ", fill = TRUE)
  }, error = function(e) {
    lut <- lookup_table
  })
  lut <- lut[,2:4]
  crp <- colorRampPalette(apply(lut, 1, function(x) {rgb(x[1], x[2], x[3], alpha)}))
  return(crp)
}
