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
               o=double(length=num.files), PACKAGE="RMINC")$o
  class(output) <- c("mincVoxel", "vector")
  attr(output, "filenames") <- filenames
  attr(output, "voxelCoord") <- c(v1,v2,v3)
  attr(output, "worldCoord") <- mincConvertVoxelToWorld(filenames[1],v1,v2,v3)
  return(output)
}

# test to see whether files exist and are readable
mincFileCheck <- function(filenames) {
  for(i in 1:length(filenames) ) {
    if(file.access(as.character(filenames[i]), 4) == -1 ){
        stop("The following file could not be read (full filename is between the dashes): ---", filenames[i], "---")
    }
  }
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
               o=double(length=num.files), PACKAGE="RMINC")$o
  class(output) <- c("mincVoxel", "vector")
  attr(output, "filenames") <- filenames
  attr(output, "worldCoord") <- c(v1,v2,v3)
  attr(output, "voxelCoord") <- mincConvertWorldToVoxel(filenames[1],v1,v2,v3)
  return(output)
}

# the print function for a voxel
print.mincVoxel <- function(x, ..., filenames=FALSE, digits=NULL) {
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
                  as.integer(v3), PACKAGE="RMINC")
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
               o=double(length=3), PACKAGE="RMINC")$o
  return(output)
}

mincConvertWorldToVoxel <- function(filename, v1, v2, v3) {
  output <- .C("convert_world_to_voxel",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3), PACKAGE="RMINC")$o
  return(round(output))
}

# return a volume as a 1D array.
mincGetVolume <- function(filename) {
  mincFileCheck(filename)
  sizes <- minc.dimensions.sizes(filename)
  start <- c(0,0,0)
  total.size <- sizes[1] * sizes[2] * sizes[3]
  output <- .C("get_hyperslab",
               as.character(filename),
               as.integer(start),
               as.integer(sizes),
               hs=double(total.size), PACKAGE="RMINC")$hs
  class(output) <- c("mincSingleDim", "numeric")
  attr(output, "filename") <- filename
  attr(output, "likeVolume") <- filename
  return(output)
}

# print function for multidimensional files
print.mincMultiDim <- function(x, ...) {
  cat("Multidimensional MINC volume\n")
  cat("Columns:      ", colnames(x), "\n")
  print(attr(x, "likeVolume"))
}

print.vertexMultiDim <- function(x, ...) {
  cat("Multidimensional Vertstats file\n")
  cat("Columns:      ", colnames(x), "\n")
}
      

print.mincSingleDim <- function(x, ...) {
  cat("MINC volume\n")
  print(attr(x, "likeVolume"))
  print(attr(x, "filename"))
}

print.mincQvals <- function(x, ...) {
  print.mincMultiDim(x)
  cat("Degrees of Freedom:", paste(attr(x, "DF")), "\n")
  cat("FDR Thresholds:\n")
  print(attr(x, "thresholds"))
}

mincWriteVolume <- function(buffer, ...) {
  UseMethod("mincWriteVolume")
}

mincWriteVolume.mincSingleDim <- function(buffer, output.filename, clobber = NULL) {
  mincWriteVolume.default(buffer, output.filename, attr(buffer, "likeVolume"), clobber)
}

# write out one column of a multidim MINC volume
mincWriteVolume.mincMultiDim <- function(buffer, output.filename, column=1, 
																				like.filename = NULL, clobber = NULL) {
  cat("Writing column", column, "to file", output.filename, "\n")
  if (is.null(like.filename)) {
    like.filename <- attr(buffer, "likeVolume")
  }
  if (is.na(file.info(as.character(like.filename))$size)) {
    stop(c("File ", like.filename, " cannot be found.\n"))
  }

  mincWriteVolume.default(buffer[,column], output.filename, like.filename, clobber)
}

# the default MINC output function
# the buffer is a vector in this case
mincWriteVolume.default <- function(buffer, output.filename, like.filename,
                                    clobber = NULL) {
	
  if(file.exists(output.filename) && is.null(clobber)){
    answer <- readline("Warning: the outputfile already exists, continue? (y/n) ")
    if(substr(answer, 1, 1) == "n")
      stop("Output file exists, specify clobber, or change the output file name.")
  }
  else if(file.exists(output.filename) && !clobber){
    stop("Output file exists, specify clobber, or change the output file name.")
  }
	
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
               as.double(buffer), PACKAGE="RMINC")
}


# get the dimension sizes of a particular file.
minc.dimensions.sizes <- function(filename) {
  sizes <- .C("get_volume_sizes",
              as.character(filename),
              sizes = integer(3), PACKAGE="RMINC")$sizes
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
               as.integer(count), hs=buffer, PACKAGE="RMINC")$hs
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

# compute a sequential ANOVA at each voxel
mincAnova <- function(formula, data=NULL, subset=NULL, mask=NULL) {
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  filenames <- as.character(mf[,1])
  mmatrix <- model.matrix(formula, mf)

  method <- "anova"

  mincFileCheck(filenames)

  v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  #l <- lm(formula, mf)
  
###   result <- .Call("minc2_model",
###                   as.character(filenames),
###                   as.matrix(mmatrix),
###                   attr(mmatrix, "assign"),
###                   as.double(! is.null(mask)),
###                   as.character(mask),
###                   NULL, NULL,
###                   as.character(method), PACKAGE="RMINC")
  result <- .Call("per_voxel_anova",
                  as.character(filenames),
                  as.matrix(mmatrix),
                  attr(mmatrix, "assign"),
                  as.double(! is.null(mask)),
                  as.character(mask), PACKAGE="RMINC")
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- rep("F", ncol(result))

  l <- lm.fit(mmatrix, v.firstVoxel)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
    dflist[[i]] <- c(dfs[[i+1]], dfr)
  }
  attr(result, "df") <- dflist

  # get the first voxel in order to get the dimension names
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))

  colnames(result) <- attr(terms(formula), "term.labels")
  class(result) <- c("mincMultiDim", "matrix")

  # run the garbage collector...
  gcout <- gc()

  return(result)
}

###########################################################################################
#' @description Linear Model at Every Voxel
#' @name mincLm
#' @title Linear model at Every Voxel
#' @param formula The linear model formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.
#' @param data The dataframe which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same length as the input files. The linear model will only be computed
#' inside the mask.
#' @details This function computes a linear model at every voxel of a set of files. The function is a close cousin to lm, the key difference
#' being that the left-hand side of the formula specification takes a series of filenames for MINC files.
#' @return mincLm returns a mincMultiDim object which contains a series of columns corresponding to the terms in the linear model. The first
#' column is the F-statistic of the significance of the entire volume, the following columns contain the marginal t-statistics for each of the terms in 
#' the model 
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples 
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(filenames ~ Genotype, gf)
#' # see what's in the results
#' vs
#' # write the results to file
#' mincWriteVolume(vs, "output.mnc", "Genotype+")
###########################################################################################
mincLm <- function(formula, data=NULL, subset=NULL, mask=NULL, maskval=NULL) {
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

  mincFileCheck(filenames)
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  result <- .Call("minc2_model",
                  as.character(filenames),
                  as.matrix(mmatrix),
                  NULL,
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.double(minmask),
                  as.double(maxmask),
                  NULL, NULL,
                  as.character(method), PACKAGE="RMINC")

  attr(result, "likeVolume") <- filenames[1]
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  
  # the order of return values is:
  #
  # f-statistic
  # r-squared
  # betas
  # t-stats
  #
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))

  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  dflist <- vector("list", (ncol(result)-2)/2 + 1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  # get the first voxel in order to get the dimension names
  v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  betaNames = paste('beta-',rows, sep='')
  tnames = paste('tvalue-',rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

# two tailed version of pt
pt2 <- function(q, df,log.p=FALSE) {
  2*pt(-abs(q), df, log.p=log.p)
}

# returns a mask as a vector - either by loading the file
# or just passing through the vector passed in.
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

vertexFDR <- function(buffer, method="FDR") {
  mincFDR.mincMultiDim(buffer, columns=NULL, mask=NULL, df=NULL,
                       method=method)
}

# mincFDR for data not created in the same R session; i.e. obtained
# from mincGetVolume
mincFDR.mincSingleDim <- function(buffer, df, mask=NULL, method="qvalue", ...) {
  if (is.null(df)) {
    stop("Error: need to specify the degrees of freedom")
  }
  if (length(df) == 1) {
    df <- c(1,df)
  }
  mincFDR.mincMultiDim(buffer, columns=1, mask=mask, df=df, method=method, ...)
}

# mincFDR for a local buffer created by mincLm
mincFDR.mincMultiDim <- function(buffer, columns=NULL, mask=NULL, df=NULL,
                                 method="FDR", statType=NULL) {
  if (method == "qvalue") {
    test <- try(library(qvalue))
    if (class(test) == "try-error") {
      stop("The qvalue package must be installed for mincFDR to work")
    }
  }


  # Remove coefficients from buffer
  stattype = attr(buffer, "stat-type")
  df  = attr(buffer,"df")
  for (nStat in 1:length(stattype)) {
    if(stattype[nStat] == 'beta' || stattype[nStat] == 'R-squared') {
      if(!exists('indicesToRemove')) {
        indicesToRemove = nStat 
      }
      else {
        indicesToRemove = c(indicesToRemove,nStat) 
      }
    }
  }
  if(exists('indicesToRemove')) {
    buffer = buffer[,-indicesToRemove]
    attr(buffer, "stat-type") <- stattype[-indicesToRemove]
  }
  attr(buffer, "df") <- df

  # must know the type of statistic we are dealing with
  knownStats <- c("t", "F")
  if (is.null(statType)) {
    # stat type not specified - must be an attribute to the buffer
    if (is.null(attr(buffer, "stat-type"))) {
      stop("Error: need to specify the type of statistic.")
    }
    else {
      statType <- attr(buffer, "stat-type")
    }
    # make sure that the stat type is recognized
    if (! all(statType %in% knownStats)) {
      stop("Error: not all the stat types are recognized. Currently allowed are: ", paste(knownStats, collapse=" "))
    }
    # make sure that there are either just one stat type
    # or as many as there are columns
    if (length(statType) == 1 & ncol(buffer) !=1) {
      statType <- rep(statType, ncol(buffer))
    }
    else if (length(statType) == ncol(buffer)) {
      # do nothing
    }
    else {
      stop("Error: stat type needs to be either a single entry or as many entries as there are columns in the buffer")
    }
  }

  # need to know the degrees of freedom
  if (is.null(df)) {
    df <- attr(buffer, "df")
    if (is.null(df)) {
      stop("Error: need to specify the degrees of freedom")
    }
    if (length(df) == 1 & ncol(buffer) != 1) {
      df <- rep(list(df), ncol(buffer))
    }
    else if (length(df) == ncol(buffer)) {
      # do nothing
    }
    else {
      stop("Error: df needs to be of either length 1 or the same length as number of columns in the buffer")
    }
    #df <- vector(length=2)
    #df[1] <- ncol(attributes(buffer)$model) -1
    #df[2] <- nrow(attributes(buffer)$model) - ncol(attributes(buffer)$model)
  }

  
  if (is.null(columns)) {
    columns <- colnames(buffer)
    cat("\nComputing FDR threshold for all columns\n")
  }

  n.cols <- length(columns)
  n.row <-0
  if (is.matrix(buffer)) {
    n.row <- nrow(buffer)
  }
  else {
    n.row <- length(buffer)
  }
  
  if (is.null(mask)) {
    mask <- vector(length=n.row) + 1
  }
  else {
    mask <- mincGetMask(mask)
  }
  

  output <- matrix(1, nrow=n.row, ncol=n.cols)
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=n.cols)
  
  for (i in 1:n.cols) {
    cat("  Computing threshold for ", columns[i], "\n")
    pvals <- 0
    qobj <- vector("list", length(pvals))

    # convert statistics to p-values
    if (statType[i] == "t") {
      if (is.matrix(buffer)) {
        pvals <- pt2(buffer[mask>0.5, i], df[[i]])
      }
      else {
        pvals <- pt2(buffer[mask>0.5], df[[i]])
      }
    }
    else if (statType[i] == "F") {
      if (is.matrix(buffer)) {
        pvals <- pf(buffer[mask>0.5, i], df[[i]][1], df[[i]][2],
                    lower.tail=FALSE)
      }
      else {
        pvals <- pf(buffer[mask>0.5], df[[i]][1], df[[i]][2], lower.tail=FALSE)
      }
    }
    



    # determine corresponding q values
    if (method=="qvalue") {
      qobj <- qvalue(pvals)
    }
    else if (method == "FDR" | method == "p.adjust") {
      qobj$pvalue <- pvals
      qobj$qvalue <- p.adjust(pvals, "fdr")
    }
    else if (method == "pFDR" | method == "fastqvalue") {
      qobj <- fast.qvalue(pvals)
    }
    # calculate thresholds at different sig levels
    for (j in 1:length(p.thresholds)) {
      if (statType[i] == "F") {
		subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
		# cat(sprintf("Number of sub-threshold F p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvalues) >= 1 ) {
			thresholds[j,i] <- qf(max(subTholdPvalues), df[[i]][1], df[[i]][2], lower.tail=FALSE)
		} else { thresholds[j,i] <- NA }
      }
      else if (statType[i] == "t") {
		subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
		#cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvalues) >= 1 ) {
			thresholds[j,i] <-qt(max(subTholdPvalues)/2, df[[i]], lower.tail=FALSE)
		} else { thresholds[j,i] <- NA }
      }
    }
    output[mask>0.5,i] <- qobj$qvalue
  }
                   
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columns
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columns
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(output)
}
   

mincMean <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="mean", maskval=maskval)
  return(result)
}

mincVar <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  return(result)
}

mincSum <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="sum", maskval=maskval)
  return(result)
}

mincSd <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  result <- sqrt(result)
  return(result)
}

#
# maskval was introduced in order to run mincSummary (and mincApply) in parralel
# another way that this argument can be used is to specify a particular label
# for which mincSummary will be used
#
mincSummary <- function(filenames, grouping=NULL, mask=NULL, method="mean", maskval=NULL) {
  mincFileCheck(filenames)
  if (is.null(grouping)) {
    grouping <- rep(1, length(filenames))
  }
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  result <- .Call("minc2_model",
                  as.character(filenames),
                  as.double(grouping)-1,
                  NULL,
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.double(minmask),
                  as.double(maxmask),
                  NULL, NULL,
                  as.character(method), PACKAGE="RMINC")
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

pMincApply <- function(filenames, function.string,
                       mask=NULL, cores=4, tinyMask=FALSE, method="snowfall") {
  # if no mask exists use the entire volume
  if (is.null(mask)) {
    maskV = mincGetVolume(filenames[1])
    nVoxels = length(maskV)
    maskV[maskV >= min(maskV)] <- as.integer(cut(seq_len(nVoxels), cores)) 
  }
  else {
    maskV <- mincGetVolume(mask)
    # optionally make the mask a fraction of the original size - for testing
    if (tinyMask!=FALSE) {
      maskV[maskV>1.5] <- 0
    }
    nVoxels <- sum(maskV>0.5)
    maskV[maskV>0.5] <- as.integer(cut(seq_len(nVoxels), cores)) 
  }
 
  maskFilename <- paste("pmincApplyTmpMask-", Sys.getpid(), ".mnc", sep="")
  mincWriteVolume(maskV, maskFilename, clobber=TRUE)
  
  pout <- list()
  
  if (method == "local") {
    stop("Lovely code ... that generates inconsistent results because something somewhere is not thread safe ...")

    library(multicore)
    library(doMC)
    library(foreach)
    registerDoMC(cores)

    # run the job spread across each core
    pout <- foreach(i=1:cores) %dopar% { mincApply(filenames, function.string,
                      mask=maskFilename, maskval=i) }
    #cat("length: ", length(pout), "\n")
  }
  else if (method == "sge") {
    library(Rsge)
    # Need to use double quotes, because both sge.submit and mincApply try to evalute the functin
    function.string = enquote(function.string)
    
    l1 <- list(length=cores)
    
    # Submit one job to the queue for each segmented brain region
    for(i in 1:cores) {
       l1[[i]]<- sge.submit(mincApply,filenames,function.string, mask=maskFilename,
                      maskval=i, packages=c("RMINC"),global.savelist=sub("\\(([A-Z]|[a-z])\\)","",function.string))
   
    }
    
    # Wait for all jobs to complete
    r1 = lapply(l1,sge.job.status)
    
    while(! all (r1 == 0)) {
      Sys.sleep(4)
      r1 = lapply(l1,sge.job.status) }
      pout <- lapply(l1, sge.list.get.result)

    function.string = eval(function.string)
  }
  else if (method == "snowfall") {
    wrapper <- function(i) {
      return(mincApply(filenames, function.string, mask=maskFilename,
                       maskval=i, reduce=FALSE))
    }
    # use all cores in the current cluster if # of cores not specified
    if (is.null(cores)) {
      cores <- length(sfSocketHosts())
    }
    
    pout <- sfLapply(1:cores, wrapper)

  }
  else {
    stop("unknown execution method")
  }

  # Need to get one voxel, x, to test number of values returned from function.string
  x <- mincGetVoxel(filenames, 0,0,0)
  test <- eval(function.string) 

  # recombine the output into a single volume
  if (length(test) > 1) {
    output <- matrix(0, nrow=length(maskV), ncol=length(test))
    class(output) <- class(pout[[1]])
    attr(output, "likeVolume") <- attr(pout[[1]], "likeVolume")
  }
  else {
    output <- maskV
  }

  for(i in 1:cores) {
    if (length(test)>1) {
      output[maskV==i,] <- pout[[i]]
    }
    else {
      currentMask = pout[[i]]
      output[maskV==i] <- currentMask[maskV==i]
    }
  }
  unlink(maskFilename)
  return(output)
}

mincApply <- function(filenames, function.string, mask=NULL, maskval=NULL, reduce=FALSE) {
 
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  # Need to get one voxel, x, to test number of values returned from function.string
  x <- mincGetVoxel(filenames, 0,0,0)
  test <- eval(function.string)
  results <- .Call("minc2_model",
                   as.character(filenames),
                   function.string,
                   NULL,
                   as.double(! is.null(mask)),
                   as.character(mask),
                   as.double(minmask),
                   as.double(maxmask),
                   .GlobalEnv,
                   as.double(length(test)),
                   as.character("eval"), PACKAGE="RMINC");

  if (length(test) > 1) {
    if (reduce==TRUE) {
      maskV <- mincGetVolume(mask)
      results <- results[maskV > (minmask-0.5) & maskV < (maxmask+0.5),]
    }
    class(results) <- c("mincMultiDim", "matrix")
  }
  else {
    if (reduce==TRUE) {
      maskV <- mincGetVolume(mask)
      results <- results[maskV > (minmask-0.5) & maskV < (maxmask+0.5)]
    }
    class(results) <- c("mincSingleDim", "numeric")
  }
  attr(results, "likeVolume") <- filenames[1]
  
  # run the garbage collector...
  gcout <- gc()
  
  return(results)
}

# use the eval interface to run mixed effect models at every vertex.
# NOTE: since it uses the eval interface it suffers from several
# important flaws:
# * can only return one statistical test
# * is numbingly, dreadfully, stupifyingly slow.

mincSlowLme <- function(filenames, fixed.effect, random.effect, mask){

  # determine the number of output variables
  x <- rnorm(length(filenames))
  s <- summary(lme(as.formula(fixed.effect), random=as.formula(random.effect)))$tTable[,4]
  l <- length(s)
  voxel.slow.lme <- function(x) {
    s <- summary(lme(as.formula(fixed.effect), random=as.formula(random.effect)))
    if (inherits(x, "try-error")) {
      return(vector("numeric", length=l))
    }
    else {
      return(s$tTable[,4])
    }
    
  }
  assign("voxel.slow.lme", voxel.slow.lme, env=.GlobalEnv)
  output <- mincApply(filenames, quote(voxel.slow.lme(x)), mask)
  return(output)

}

mincLme <- function(filenames, fixed.effect, random.effect, mask=NULL)
{
  # determine the number of output variables
  x <- rnorm(length(filenames))
  s <- rmincLme(as.formula(fixed.effect), random=as.formula(random.effect))

  voxel.lme <- function(x) {
    s <- rmincLmeLoop(dataMix, X, Z, grps, lmeSt, controlvals, dims, listNncols, listrownames, x)
    if (inherits(s, "try-error")) {
      return(vector("numeric", length=2))
    }
    else {
	return(s)
    }
    
  }
  assign("voxel.lme", voxel.lme, env=.GlobalEnv)
  output <- mincApplyLme(filenames, quote(voxel.lme(x)), mask)
  return(output)

}

mincApplyLme <- function(filenames, function.string, mask=NULL, maskval=NULL) {
  x <- mincGetVoxel(filenames, 0,0,0)
  test <- eval(function.string)
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  results <- .Call("minc2_model",
                   as.character(filenames),
                   function.string,
                   NULL,
                   as.double(! is.null(mask)),
                   as.character(mask),
                   as.double(minmask),
                   as.double(maxmask),
                   .GlobalEnv,
                   as.double(length(test)),
                   as.character("eval"), PACKAGE="RMINC");

  attr(results, "likeVolume") <- filenames[1]
  if (length(test) > 1) {
    class(results) <- c("mincMultiDim", "matrix")
  }
  else {
    class(results) <- c("mincSingleDim", "numeric")
  }
  colnames(results) <- colnames(test)
  
  # run the garbage collector...
  gcout <- gc()
  
  return(results)
}
  
vertexTable <- function(filenames) {
  nSubjects <- length(filenames)
  nvertices <- nrow(read.table(filenames[1]))
  output <- matrix(nrow=nvertices, ncol=nSubjects)
  for (i in 1:nSubjects) {
    output[,i] <- as.matrix(read.table(filenames[i]))
  }
  return(output)
}

###########################################################################################
#' Performs ANOVA on each vertex point specified 
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param filenames list of vertex files
#' @param subset rows to be used, by default all are used
#' @return Returns an array with the F-statistic for each model specified by formula with the following attributes: model – design matrix, filenames – 
#' 	vertex file names input, stat-type: type of statistic used, df – degrees of freedom of each statistic. 
#' @seealso mincAnova,anatAnova 
#' @examples 
#' gf = read.csv("~/SubjectTable.csv") 
#' civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
#' gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
#' result = vertexAnova(~Primary.Diagnosis,gf,gf$CIVETFILES$nativeRMStlink20mmleft) 
###########################################################################################
vertexAnova <- function(formula, data=NULL,filenames, subset=NULL) {
  # Create Model
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mmatrix <- model.matrix(formula, mf)

  # Load Vertex Data from Files
  #filenames <- as.character(mf[,1])
  data.matrix <- vertexTable(filenames)
  result <- .Call("vertex_anova_loop", data.matrix, mmatrix,attr(mmatrix, "assign"), PACKAGE="RMINC");

  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- rep("F", ncol(result))

  # Use the first voxel in order to get the dimension names
  v.firstVoxel <- data.matrix[1,]
  columns <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  assignVector = attr(mmatrix, "assign") + 1
  columnName =  rep('', max(assignVector)-1)
  dflist =  rep(0, max(assignVector)-1)
  for (i in 2:max(assignVector)) { 
    indices = which(assignVector == i)
    for (j in 1:length(indices)) {
      columnName[i-1] = paste(columnName[i-1],'',columns[indices[j]]) 
    }
    dflist[i-1] = length(indices)
  }

  attr(result, "df") <- dflist
  colnames(result) <- columnName

  class(result) <- c("vertexMultiDim", "matrix")
  return(result)
}
###########################################################################################
#' Calculates statistics and coefficients for linear model of specified vertex files
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param subset rows to be used, by default all are used
#' @return Returns an object containing the beta coefficients, F 
#' and t statistcs that can be passed directly into vertexFDR.
#' @seealso mincLm,anatLm,vertexFDR 
#' @examples 
#' gf = read.csv("~/SubjectTable.csv") 
#' civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
#' gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
#' gf$vertexFiles = as.factor(gf$CIVETFILES$nativeRMStlink20mmleft)
#' result = vertexLm(vertexFiles~Primary.Diagnosis,gf) 
#' vertexFDR(result)
###########################################################################################
vertexLm <- function(formula, data, subset=NULL) {
  # repeat code to extract the formula as in mincLm
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)

  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())


  filenames <- as.character(mf[,1])
  mmatrix <- model.matrix(formula, mf)

  cat("Loading data from files\n")
  data.matrix <- vertexTable(filenames)

  cat("after loading\n")

  result <- .Call("vertex_lm_loop", data.matrix, mmatrix, PACKAGE="RMINC");

  attr(result, "likeVolume") <- filenames[1]
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
 
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2+1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  # get the first voxel in order to get the dimension names
  v.firstVoxel <- data.matrix[1,]
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))

  betaNames = paste('beta-', rows, sep='')
  tnames = paste('tvalue-', rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("vertexMultiDim", "matrix")
 
  # run the garbage collector...
  gcout <- gc()
 
  return(result)
}
###########################################################################################
#' Writes vertex data to a file with an optional header
#' @param vertexData vertex data to be written
#' @param filename full path to file where data shall be written
#' @param headers Whether or not to write header information
#' @param mean.stats mean vertex data that may also be written
#' @param gf glim matrix that can be written to the header
#' @return A file is generated with the vertex data and optional headers
#' @examples 
#' gf = read.csv("~/SubjectTable.csv") 
#' gfCIVET = civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
#' gfCIVET = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gfCIVET)
#' writeVertex(gfCIVET$nativeRMStlink20mm,"~/RMStlink20mm",TRUE,NULL,gf)
###########################################################################################
writeVertex <- function (vertexData, filename, headers = TRUE, mean.stats = NULL, 
    gf = NULL) 
{
    append.file = TRUE
    if (headers == TRUE) {
        write("<header>", file = filename)
        if (is.object(mean.stats)) {
            write("<mean>", file = filename, append = TRUE)
            sink(filename, append = TRUE)
            print(summary(mean.stats))
            sink(NULL)
            write("</mean>", file = filename, append = TRUE)
            write("<formula>", file = filename, append = TRUE)
            sink(filename, append = TRUE)
            print(formula(mean.stats))
            sink(NULL)
            write("</formula>", file = filename, append = TRUE)
        }
        if (is.data.frame(glim.matrix)) {
            write("<matrix>", file = filename, append = TRUE)
            write.table(glim.matrix, file = filename, append = TRUE, 
                row.names = FALSE)
            write("</matrix>", file = filename, append = TRUE)
        }
        write("</header>", file = filename, append = TRUE)
    }
    else {
        append.file = FALSE
    }
    write.table(vertexData, file = filename, append = append.file, 
        quote = FALSE, row.names = FALSE, col.names = headers)
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

# calls ray_trace_crosshair.pl to generate a pretty picture of an anatomical slice
# with a stats slice overlayed and a crosshair indicating the coordinate (peak)
# of interest
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
