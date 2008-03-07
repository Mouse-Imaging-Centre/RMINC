# get the real value of one voxel from all files.
mincGetVoxel <- function(filenames, v1, v2=NULL, v3=NULL) {
  num.files <- length(filenames)
  if (length(v1) == 3){
    v2 <- v1[2]
    v3 <- v1[3]
    v1 <- v1[1]
  }
  else if (is.null(v2) || is.null(v3)) {
    stop("Three elements have to be specified.")
  }
  output <- .C("get_voxel_from_files",
               as.character(filenames),
               as.integer(num.files),
               as.integer(v1),
               as.integer(v2),
               as.integer(v3),
               o=double(length=num.files))$o
  class(output) <- c("mincVoxel", class(output))
  attr(output, "filenames") <- filenames
  attr(output, "voxelCoord") <- c(v1,v2,v3)
  attr(output, "worldCoord") <- mincConvertVoxelToWorld(filenames[1],v1,v2,v3)
  return(output)
}

# get the real value of one voxel from all files using world coordinates
mincGetWorldVoxel <- function(filenames, v1, v2=NULL, v3=NULL) {
  num.files <- length(filenames)
  if (length(v1) == 3){
    v2 <- v1[2]
    v3 <- v1[3]
    v1 <- v1[1]
  }
  else if (is.null(v2) || is.null(v3)) {
    stop("Three elements have to be specified.")
  }
  output <- .C("get_world_voxel_from_files",
               as.character(filenames),
               as.integer(num.files),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=num.files))$o
  class(output) <- c("mincVoxel", class(output))
  attr(output, "filenames") <- filenames
  attr(output, "worldCoord") <- c(v1,v2,v3)
  attr(output, "voxelCoord") <- mincConvertWorldToVoxel(filenames[1],v1,v2,v3)
  return(output)
}

# the print function for a voxel
print.mincVoxel <- function(x, filenames=FALSE) {
  if (filenames == FALSE) {
    print.table(x)
  }
  else {
    print.table(cbind(x, attr(x,"filenames")))
  }
  cat("\nVoxel Coordinates:", attr(x, "voxelCoord"), "\n")
  cat("World Coordinates:", attr(x, "worldCoord"), "\n")
}

# gets a vector from a series of 4D minc volumes
mincGetVector <- function(filenames, v1, v2, v3, v.length) {
  num.files <- length(filenames)
  output <- .Call("get_vector_from_files",
                  as.character(filenames),
                  as.integer(num.files),
                  as.integer(v.length),
                  as.integer(v1),
                  as.integer(v2),
                  as.integer(v3))
  class(output) <- c("mincVector", "mincVoxel", class(output))
  attr(output, "filenames") <- filenames
  attr(output, "voxelCoord") <- c(v1,v2,v3)
  attr(output, "worldCoord") <- mincConvertVoxelToWorld(filenames[1],v1,v2,v3)
  return(output)
}

mincConvertVoxelToWorld <- function(filename, v1, v2, v3) {
  output <- .C("convert_voxel_to_world",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3))$o
  return(output)
}

mincConvertWorldToVoxel <- function(filename, v1, v2, v3) {
  output <- .C("convert_world_to_voxel",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3))$o
  return(round(output))
}

# return a volume as a 1D array.
mincGetVolume <- function(filename) {
  sizes <- minc.dimensions.sizes(filename)
  start <- c(0,0,0)
  total.size <- sizes[1] * sizes[2] * sizes[3]
  output <- .C("get_hyperslab",
               as.character(filename),
               as.integer(start),
               as.integer(sizes),
               hs=double(total.size))$hs
  class(output) <- c("mincSingleDim", "numeric")
  return(output)
}

# print function for multidimensional files
print.mincMultiDim <- function(x) {
  cat("Multidimensional MINC volume\n")
  cat("Columns:      ", colnames(x), "\n")
  print(attr(x, "likeVolume"))
}

print.mincSingleDim <- function(x) {
  cat("MINC volume\n")
  print(attr(x, "likeVolume"))
}

print.mincQvals <- function(x) {
  print.mincMultiDim(x)
  cat("Degrees of Freedom:", attr(x, "DF"), "\n")
  cat("FDR Thresholds:\n")
  print(attr(x, "thresholds"))
}

mincWriteVolume <- function(buffer, ...) {
  UseMethod("mincWriteVolume")
}

mincWriteVolume.mincSingleDim <- function(buffer, output.filename) {
  mincWriteVolume.mincMultiDim(buffer, output.filename)
}

# write out one column of a multidim MINC volume
mincWriteVolume.mincMultiDim <- function(buffer, output.filename, column=1, like.filename = NULL) {
  cat("Writing column", column, "to file", output.filename, "\n")
  if (is.null(like.filename)) {
    like.filename <- attr(buffer, "likeVolume")
  }
  if (is.na(file.info(like.filename)$size)) {
    stop(c("File ", like.filename, " cannot be found.\n"))
  }

  mincWriteVolume.default(buffer[,column], output.filename, like.filename)
}

# the default MINC output function
# the buffer is a vector in this case
mincWriteVolume.default <- function(buffer, output.filename, like.filename) {
  sizes <- minc.dimensions.sizes(like.filename)
  start <- c(0,0,0)
  if ((sizes[1] * sizes[2] * sizes[3]) != length(buffer)) {
    stop("Size of like-file not the same as size of buffer")
  }
  b.min <- min(buffer)
  b.max <- max(buffer)
  output <- .C("write_minc2_volume",
               as.character(output.filename),
               as.character(like.filename),
               as.integer(start),
               as.integer(sizes),
               as.double(b.max),
               as.double(b.min),
               as.double(buffer))
}


# get the dimension sizes of a particular file.
minc.dimensions.sizes <- function(filename) {
  sizes <- .C("get_volume_sizes",
              as.character(filename),
              sizes = integer(3))$sizes
  return(sizes)
}

# get a hypeslab from an existing volume with given start and counts.
minc.get.hyperslab <- function(filename, start, count, buffer=NA) {
#  total.size <- (count[1] - start[1]) * (count[2] - start[2]) *
#                (count[3] - start[3])
  total.size <- (count[1] * count[2] * count[3])
  if (is.na(buffer)) {
    buffer <- double(total.size)
  }
  cat("Total size: "); cat(total.size); cat("\n")
  output <- .C("get_hyperslab",
               as.character(filename),
               as.integer(start),
               as.integer(count), hs=buffer)$hs
  return(output)
}

# some rather MICe specific code best ignored.
wilcox.permutation.full <- function(filenames, groupings, mask, n.permute=10) {
  results <- matrix(nrow=n.permute, ncol=55)
  mask.volume <- minc.get.volume(mask)
  for (i in 1:n.permute) {
    new.order <- sample(groupings)
    cat("N: ", as.double(new.order), "\n")
    w <- minc.wilcoxon.test(filenames, new.order, mask)
    w2 <- w[mask.volume == 1]
    results[i,] <- tabulate(round(w2),55)
  }
  return(results)
}

f <- function(formula, data=NULL, subset=NULL, mask=NULL) {
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  filenames <- mf[,1]
  mmatrix <- model.matrix(formula, mf)
  
  cat("MASK: ", mask, "\n")
  
  return(mmatrix)
}



mincLm <- function(formula, data=NULL, subset=NULL, mask=NULL) {
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  filenames <- as.character(mf[,1])
  mmatrix <- model.matrix(formula, mf)

  method <- "lm"

  result <- .Call("minc2_model",
                  as.character(filenames),
                  as.matrix(mmatrix),
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.character(method))
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames

  # get the first voxel in order to get the dimension names
  v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))

  colnames(result) <- c("F-statistic", rows)
  class(result) <- c("mincMultiDim", "matrix")

  return(result)
}

pt2 <- function(q, df,log.p=FALSE) {
  2*pt(-abs(q), df, log.p=log.p)
}

mincGetMask <- function(mask) {
  if (class(mask) == "character") {
    return(mincGetVolume(mask))
  }
  else {
    return(mask)
  }
}

mincFDR <- function(buffer, ...) {
  UseMethod("mincFDR")
}

mincFDR.mincMultiDim <- function(buffer, columns=NULL, mask=NULL, df=NULL) {
  test <- try(library(qvalue))
  if (class(test) == "try-error") {
    stop("The qvalue package must be installed for mincFDR to work")
  }

  if (is.null(mask)) {
    mask <- vector(length=nrow(buffer)) + 1
  }
  else {
    mask <- mincGetMask(mask)
  }
  
  if (is.null(df)) {
    df <- vector(length=2)
    df[1] <- ncol(attributes(buffer)$model) -1
    df[2] <- nrow(attributes(buffer)$model) - ncol(attributes(buffer)$model)
  }

  if (is.null(columns)) {
    columns <- colnames(buffer)
    cat("\nComputing FDR threshold for all columns\n")
  }

  n.cols <- length(columns)
  output <- matrix(1, nrow=nrow(buffer), ncol=n.cols)
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=n.cols)
  
  for (i in 1:n.cols) {
    cat("  Computing threshold for ", columns[i], "\n")
    qobj <- 0
    if (columns[i] == "F-statistic") {
      qobj <- qvalue(pf(buffer[mask > 0.5, columns[i]],
                        df[1], df[2], lower.tail=F))
      for (j in 1:length(p.thresholds)) {
        thresholds[j,i] <- qf(max(qobj$pvalue[qobj$qvalue <= p.thresholds[j]]),
                                  df[1], df[2], lower.tail=F)
      }
      output[mask>0.5,i] <- qobj$qvalue
    }
                   
    else {
      qobj <- qvalue(pt2(buffer[mask>0.5, columns[i]], df[2]))
      for (j in 1:length(p.thresholds)) {
        thresholds[j,i] <-qt(max(qobj$pvalue[qobj$qvalue <= p.thresholds[j]])/2,
                              df[2], lower.tail=F)
      }
      output[mask>0.5,i] <- qobj$qvalue
    }
  }
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columns
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columns
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim")
  return(output)
}
   

mincMean <- function(filenames, grouping=NULL, mask=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="mean")
  return(result)
}

mincVar <- function(filenames, grouping=NULL, mask=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var")
  return(result)
}

mincSum <- function(filenames, grouping=NULL, mask=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="sum")
  return(result)
}

mincSd <- function(filenames, grouping=NULL, mask=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var")
  result <- sqrt(result)
  return(result)
}


mincSummary <- function(filenames, grouping=NULL, mask=NULL, method="mean") {
  if (is.null(grouping)) {
    grouping <- rep(1, length(filenames))
  }

  result <- .Call("minc2_model",
                  as.character(filenames),
                  as.double(grouping)-1,
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.character(method))
  attr(result, "likeVolume") <- as.character(filenames[1])
  attr(result, "filenames") <- as.character(filenames)


  if (is.null(grouping)) {
    class(result) <- c("mincSingleDim", "numeric")
  }
  else {
    class(result) <- c("mincMultiDim", "matrix")
    colnames(result) <- levels(grouping)
  }
  return(result)
}
  
# run a t-test, wilcoxon test, correlation, or linear model at every voxel
minc.model <- function(filenames, groupings, method="t-test",
                       mask=NULL) {

  if (method == "t-test"
      || method == "wilcoxon"
      || method == "correlation"
      || method == "lm"
      || method == "paired-t-test") {
    # do nothing
  }
  else {
    stop("Method must be one of t-test, paired-t-test, wilcoxon, correlation or lm")
  }

  if (method == "lm") {
    result <- list(method="lm")
    result$likeVolume <- filenames[1]
    result$model <- groupings
    result$filenames <- filenames
    result$data <- .Call("minc2_model",
                         as.character(filenames),
                         as.matrix(groupings),
                         as.double(! is.null(mask)),
                         as.character(mask),
                         as.character(method))

    # get the first voxel in order to get the dimension names
    v.firstVoxel <- minc.get.voxel.from.files(filenames, 0,0,0)
    rows <- sub('mmatrix', '',
                rownames(summary(lm(v.firstVoxel ~ groupings))$coefficients))
    colnames(result$data) <- c("F-statistic", rows)
    class(result) <- "mincMultiDim"
    
  }
  else {
    groupings <- as.double(groupings)
    result <- .Call("minc2_model",
                    as.character(filenames),
                    as.double(groupings),
                    as.double(! is.null(mask)),
                    as.character(mask),
                    as.character(method))
  }
  return(result)
}

# create a 2D array of full volumes of all files specified.
minc.get.volumes <- function(filenames) {
  sizes <- minc.dimensions.sizes(filenames[1])
  n.files <- length(filenames)
  output <- matrix(ncol=n.files, nrow=(sizes[1] * sizes[2] * sizes[3]))
  for (i in 1:n.files) {
    output[,i] <- minc.get.volume(filenames[i])
  }
  return(output)
}


# efficient way of applying an R function to every voxel
mincApply <- function(filenames, function.string, mask=NULL) {
  result <- list(method=paste("mincApply:", function.string))
  result$likeVolume <- filenames[1]
  result$filenames <- filenames
  result$data <- .Call("minc2_apply", as.character(filenames),
                       function.string,
                       as.double(! is.null(mask)),
                       as.character(mask),
                       parent.env(environment()))
}

# calls ray-trace to generate a pretty picture of a slice
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
    vol <- minc.get.volume(volume)
    threshold <- range(vol)
    rm(vol)
  }

  # get the background threshold if necessary
  if (!is.null(background) && is.null(background.threshold)) {
    vol <- minc.get.volume(background)
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


# use the eval interface to run mixed effect models at every vertex.
# NOTE: since it uses the eval interface it suffers from several
# important flaws:
# * can only return one statistical test
# * is numbingly, dreadfully, stupifyingly slow.

minc.slow.lme <- function(filenames, fixed.effect, random.effect,
                          column, mask){
  voxel.slow.lme <- function(x) {
    summary(lme(as.formula(fixed.effect), random=as.formula(random.effect)))$tTable[column,4]
  }
  assign("voxel.slow.lme", voxel.slow.lme, env=.GlobalEnv)
  output <- minc.apply(filenames, quote(voxel.slow.lme(x)), mask)
  return(output)

}
  

# calls ray_trace_crosshair.pl to generate a pretty picture of an anatomical slice
# with a stats slice overlayed and a crosshair indicating the coordinate (peak)
# of interest
mincRayTraceStats <- function(v, anatomy.volume, stats, caption="t-statistic",
fdr=NULL, slice.direction="transverse",
outputfile="ray_trace_crosshair.png", show.pos.and.neg=FALSE, display=TRUE) 
{
	#check whether ray_trace_crosshair is installed
	lasterr <- try(system("ray_trace_crosshair", ignore.stderr = TRUE), silent=TRUE)
	if(lasterr == 32512)
	{
		stop("ray_trace_crosshair must be installed for mincRayTraceStats to work.")
	}
	
	#check whether the first argument is a mincVoxel:
	if(class(v)[1] != "mincVoxel")
	{
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
	
	i <- 10;
	if(class(anatomy.volume)[1] == "character")
	{
		systemcall[i] <- "-final-nlin"
		i <- i + 1
		systemcall[i] <- anatomy.volume
		i <- i + 1
	}
	else #(class(anatomy.volume)[1] == "mincSingleDim")
	{
		stop("Please specify the path to the anatomy volume.")
	}
	
	if(class(stats) == "character")
	{
		systemcall[i] <- "-jacobian"
		i <- i + 1
		systemcall[i] <- stats
		i <- i + 1
	}
	else if(class(stats) == "numeric")
	{
		#write buffer to file
		mincWriteVolume(stats, "/tmp/R-wrapper-ray-trace-stats.mnc", anatomy.volume)
		systemcall[i] <- "-jacobian"
		i <- i + 1
		systemcall[i] <- "/tmp/R-wrapper-ray-trace-stats.mnc"
		i <- i + 1
	}
	
	if(class(fdr) == "numeric")
	{
		systemcall[i] <- "-fdr"
		i <- i + 1
		systemcall[i] <- fdr
		i <- i + 1
	}
	
	systemcall[i] <- "-outputfile"
	i <- i + 1
	systemcall[i] <- outputfile
	i <- i + 1
	
	systemcall[i] <- "-slicedirection"
	i <- i + 1
	systemcall[i] <- slice.direction
	i <- i + 1
	
	if(show.pos.and.neg == FALSE)
	{
		systemcall[i] <- "-no-show-pos-and-neg"
		i <- i + 1
	}
	else
	{
		systemcall[i] <- "-show-pos-and-neg"
		i <- i + 1
	}

	systemcall[i] <- "-remove-temp"
	i <- i + 1
	
 	system(paste(systemcall, collapse = " " ))
 						
  if(display) 
  {
    system(paste("display", outputfile, "&"))
  }
  
  if(class(stats) ==  "numeric")
	{
		#remove the file written to disk
		system(paste("rm -f /tmp/R-wrapper-ray-trace-stats.mnc"))
	}
  
}

