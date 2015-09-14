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
                                begin=NULL, end=NULL, symmetric=F) {
  opar <- par()
  par(mfrow=mfrow)
  nslices <- prod(mfrow)
  par(mar=c(0,0,0,0))
  d <- dim(anatomy)
  if (is.null(begin)) { begin <- 1 }
  if (is.null(end)) { end <- d[dimension]-1 }
  else if (end < 0) { end <- d[dimension] + end }
  slices <- ceiling(seq(begin, end, length=nslices))
  if (is.null(anatLow)) { anatLow <- hist(anatomy, plot=F)$mid[5] }
  if (is.null(anatHigh)) { anatHigh <- rev(hist(anatomy, plot=F)$mid)[5]}
  for (i in 1:nslices) {
    mincPlotAnatAndStatsSlice(anatomy, statistics, dimension, slice=slices[i],
                              low=low, high=high, anatLow=anatLow, anatHigh=anatHigh,
                              col=col, legend=NULL, symmetric=symmetric)
  }
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
  a <- getSlice(anatomy, slice, dimension)
  s <- getSlice(statistics, slice, dimension)
  image(scaleSlice(a$slice, anatLow, anatHigh, underTransparent=F), 
        useRaster=T, col=gray.colors(100), asp=a$asp, axes=F)
  image(scaleSlice(s$slice, low, high), useRaster=T, col=col, asp=s$asp, 
        add=T, axes=F)
  if (symmetric==TRUE) {
    image(scaleSlice(s$slice, low*-1, high*-1), useRaster=T, col=rcol, asp=s$asp, 
          add=T, axes=F)
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