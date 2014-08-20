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

# 
###########################################################################################
#' @description Load a 3-dimensional MINC2 volume and returns it as a 1D array
#' @name mincGetVolume
#' @title Return a volume as a 1D array
#' @param filename A string corresponding to the location of the MINC2 file to
#' be read.
#' @return mincLm Returns a vector of mincSingleDim class
#' @seealso mincWriteVolume
###########################################################################################
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

print.mincLogLikRatio <- function(x, ...) {
  cat("mincLogLikRatio output\n\n")
  mincLmerLists <- attr(x, "mincLmerLists")
  for (i in 1:length(mincLmerLists)) {
    cat("Model", i, ":\n")
    cat("  Formula:  ")
    print(mincLmerLists[[i]][[1]]$formula)
  }
  cat("\nMask used:", attr(out, "mask"), "\n")
  cat("Chi-squared Degrees of Freedom:", attr(x, "df"), "\n")
}

print.mincLmer <- function(x, ...) {
  cat("mincLmer output\n")
  cat("Formula:  ")
  print(attr(x, "mincLmerList")[[1]]$formula)
  if (attr(x,"mincLmerList")[[1]]$REML == TRUE) {
    cat("Fitted with REML\n")
  }
  else {
    cat("Fitted with ML\n")
  }
  cat("Mask:         ", attr(x, "mask"), "\n")
  cat("Columns:      ", colnames(x), "\n")
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


###########################################################################################
#' @description Writes a MINC volume to file
#' @name mincWriteVolume
#' @aliases mincWriteVolume.mincSingleDim,mincWriteVolume.mincMultiDim
#' @title Write a MINC volume to file
#' @param buffer The data to be written to file. Usually the result of
#' \link{mincLm} or some such command
#' @param output.filename The filename to which to write the data to
#' @param column Optional name of the column of a multidimensional MINC
#' object to write out. By default the first column is used
#' @param like.filename An existing MINC filename which has the same
#' dimensions as the data to be written out. Normally this information
#' is stored inside MINC data objects
#' @param clobber Overwrite existing output file when set to TRUE, will not 
#' overwrite when set to FALSE and will prompt when NULL
#' @details This function takes numeric data, usually the results computed
#' from one of the other mincFunctions, and writes it to file so that it
#' can be viewed or manipulated with the standard MINC tools
#' @return mincLm Returns a vector of mincSingleDim class
#' @seealso mincWriteVolume,mincLm,mincFDR,mincMean,mincSd
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
###########################################################################################

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

###########################################################################################
#' @description compute a sequential ANOVA at each voxel
#' @name mincAnova
#' @title Anova at Every Voxel
#' @param formula The anova formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.
#' @param data The dataframe which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same length as the input files. ANOVA will only be computed
#' inside the mask.
#' @details This function computes a sequential ANOVA over a set of files.
#' @return Returns an array with the F-statistic for each model specified by formula with the following attributes: model – design matrix, filenames – 
#' 	vertex file names input, stat-type: type of statistic used, df – degrees of freedom of each statistic. 
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples 
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # run aan ANOVA relating the data in all voxels to Genotype
#' vs <- mincLm(filenames ~ Genotype, gf)
#' # see what's in the results
#' vs
###########################################################################################

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
#' @param formula The linear model formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.The RHS of the formula may contain one term with filenames. If so only the + operator may be used, and only two terms may appear on the RHS
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
mincLm <- function(formula, data=NULL,subset=NULL , mask=NULL, maskval=NULL) {

  #INITIALIZATION
  method <- "lm"
   
  # Build model.frame
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())


  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }


  # the following returns:
  #
  # list(data.matrix.left = data.matrix.left, 
  #      data.matrix.right = data.matrix.right,
  #      rows = rows,
  #      matrixFound = matrixFound,
  #      mmatrix = mmatrix)
  parseLmOutput <- parseLmFormula(formula,data,mf)

  # Call subroutine based on whether matrix was found
  if(parseLmOutput$matrixFound) {
    mincFileCheck(parseLmOutput$data.matrix.left)
    mincFileCheck(parseLmOutput$data.matrix.right)
  }
  else  {
    parseLmOutput$mmatrix <- model.matrix(formula, mf)	
    parseLmOutput$data.matrix.left <- as.character(mf[,1])
    mincFileCheck(parseLmOutput$data.matrix.left)
    parseLmOutput$rows = colnames(parseLmOutput$mmatrix)
  }

  result <- .Call("minc2_model",
                  as.character(parseLmOutput$data.matrix.left),
                  parseLmOutput$data.matrix.right,
                  as.matrix(parseLmOutput$mmatrix),
                  NULL,
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.double(minmask),
                  as.double(maxmask),
                  NULL, NULL,
                  as.character(method), PACKAGE="RMINC")

  attr(result, "likeVolume") <- parseLmOutput$data.matrix.left[1]
  attr(result, "filenames") <- parseLmOutput$data.matrix.left
  attr(result, "model") <- as.matrix(parseLmOutput$mmatrix)
 

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
  #v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  #rows <- sub('mmatrix', '',
  #            rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  betaNames = paste('beta-',parseLmOutput$rows, sep='')
  tnames = paste('tvalue-',parseLmOutput$rows, sep='')
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

###########################################################################################
#' @description Takes the output of a mincLm run and computes the False Discovery Rate on the results.
#' @name mincFDR
#' @aliases mincFDR,mincFDR.mincSingleDim,mincFDR.mincMultiDim,vertexFDR,anatFDR
#' @title Compute the False Discovery Rate for a mincLm object
#' @usage \method{mincFDR}{mincSingleDim}(buffer, df, mask=NULL, method="qvalue", \dots)
#' 	  \method{mincFDR}{mincMultiDim}(buffer, columns=NULL, mask=NULL, df=NULL,method="FDR", statType=NULL)
#' @param buffer The results of a mincLm run.
#' @param columns A vector of column names. By default the threshold will
#' be computed for all columns; with this argument the computation can
#' be limited to a subset.
#' @param mask Either a filename or a numeric vector representing a mask
#' only values inside the mask will be used to compute the
#' threshold.
#' @param df The degrees of freedom - normally this can be determined
#' from the input object.
#' @param method The method used to compute the false discovery
#' rate. Options are "FDR" and "pFDR".
#' @param statType This should be either a "t" or an "F", depending upon the type
#' of statistic being thresholded.
#' @details This function uses the \code{qvalue} package to compute the
#'  False Discovery Rate threshold for the results of a \link{mincLm}
#'  computation. The False Discovery Rate represents the percentage of
#'  results expected to be a false positive. Two implementations can be
#'  used as specified by the method argument. "FDR" uses the
#'  implementation in \code{p.adjust}, whereas "pFDR" is a version of the
#'  postivie False Discovery Rate as found in John Storey's \code{qvalue}
#'  package. 
#' @return qvals mincFDR returns an object with the same number of columns
#'  as the input (or the subset specified by the columns argument to
#'  mincFDR). Each column now contains the qvalues for each voxel. Areas
#'  outside the mask (if a mask was specified) will be represented by a
#'  value of 1. The result also has an attribute called "thresholds"
#'  which contains the 1, 5, 10, 15, and 20 percent false discovery rate
#'  thresholds.
#' @seealso mincWriteVolume,mincLm
#' @examples 
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(filenames ~ Genotype, gf)
#' # compute the False Discovery Rate
#' qvals <- mincFDR(vs, mask="mask.mnc")
#' # write the Gentoype column of the qvals to file
#' mincWriteVolume(qvals, "FDR-results.mnc", "Genotype+")
###########################################################################################
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


mincFDR.mincLmer <- function(buffer, columns=NULL, mask=NULL) {
  cat("In mincFDR.mincLmer\n")

  # if no DF set, exit with message
  df <- attr(buffer, "df")
  if (is.null(df)) {
    stop("No degrees of freedom for mincLmer object. Needs to be explicitly assigned with mincLmerEstimateDF (and read the documentation of that function to learn about the dragons that be living there!).")
  }
  else {
    warning("Here be dragons! Hypothesis testing with mixed effects models is challenging, since nobody quite knows how to correctly estimate denominator degrees of freedom.")
  }

  # NOTE: all this mask stuff could likely be encapsulated in its own function
  # check if mask was passed in, if not, use mask from buffer attribute
  if (is.null(mask)) {
    mask <- attr(buffer, "mask")
  }
  cat("Using mask:", mask, "\n")

  # if mask is still null, create a vector of ones with length of buffer
  if (is.null(mask)) {
    mask <- vector(length=nrow(buffer)) + 1
  }
  else {
    mask <- mincGetMask(mask)
  }
  
  # only compute stats on tlmer columns
  tlmerColumns <- grep("tlmer", attr(buffer, "stat-type"))
  ncolsToUse <- length(tlmerColumns)
  # sanity check to ensure that number of tlmer columns matches DF
  if (ncolsToUse != length(df)) {
    stop("Mismatch between DF and number of columns")
  }

  # compute p and q values
  # pvals and qvals size corresponds to voxels inside mask, whereas output
  # is the size of the buffer.
  pvals <- matrix(nrow=sum(mask>0.5), ncol=ncolsToUse)
  qvals <- pvals
  output <- matrix(1, nrow=nrow(buffer), ncol=ncolsToUse)

  # compute the thresholds at several sig levels
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=ncolsToUse)
  for (i in 1:ncolsToUse) {
    # compute qvals through R's p.adjust function
    pvals[, i] <- pt2(buffer[mask>0.5, tlmerColumns[i]], df[[i]])
    qvals[, i] <- p.adjust(pvals[, i], "fdr")
    output[mask>0.5, i] <- qvals[, i]
    for (j in 1:length(p.thresholds)) {
      # compute thresholds; to be honest, not quite sure what the NA checking is about
      subTholdPvalues <- pvals[qvals[,i] <= p.thresholds[j], i]
      subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
                                        
      if ( length(subTholdPvaluesNumbers) >= 1 ) {
        thresholds[j,i] <-qt(max(subTholdPvaluesNumbers)/2, df[[i]], lower.tail=FALSE)
      }
      else { thresholds[j,i] <- NA }
    }
  }

  columnNames <- colnames(buffer)[tlmerColumns]
  columnNamesQ <- sub("tvalue", "qvalue", columnNames)
  
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columnNames
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columnNamesQ
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()

  return(output)
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
    if(stattype[nStat] == 'beta' || stattype[nStat] == 'R-squared' || stattype[nStat] == "logLik") {
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
    attr(buffer, "df") <- df
  }


  # must know the type of statistic we are dealing with
  knownStats <- c("t", "F", "u", "chisq", "tlmer")
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
      stop("Error: not all the stat types are recognized. Currently allowed are: ",
           paste(knownStats, collapse=" "))
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


  if ( any(statType %in% "u")) {
    m <- attr(buffer, "m") 
    n <- attr(buffer, "n")
  }
  else {
    # need to know the degrees of freedom
    df <- attr(buffer, "df")
    if (is.null(df)) {
      if (any(statType %in% "tlmer")) {
        stop("Error: no degrees of freedom for mincLmer object. Needs to be explicitly assigned with mincLmerEstimateDF (and read the documentation of that function to learn about the dragons that be living there!).")
      }
      else {
        stop("Error: need to specify the degrees of freedom")
      }
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

  if (any("tlmer" %in% statType)) {
    warning("Warning: computing p-values from a mincLmer call. Mixed-effects models are notoriously difficult to correctly obtain p-values from, so this is based on an approximation and might be incorrect. Read the documentation and, if in doubt, use log likelihood testing for a more correct approach.")
  }
  
  for (i in 1:n.cols) {
    cat("  Computing threshold for ", columns[i], "\n")
    pvals <- 0
    qobj <- vector("list", length(pvals))

    # convert statistics to p-values
    if (statType[i] %in% c("t", "tlmer")) {
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
    else if (statType[i] == "u") {
      pvals <- 1 - pwilcox(buffer[mask>0.5,i],m,n,lower.tail = FALSE)
    }
    else if (statType[i] == "chisq") {
      if (is.matrix(buffer)) {
        pvals <- pchisq(buffer[mask>0.5, i], df[[i]], lower.tail=F)
      }
      else {
        pvals <- pchisq(buffer[mask>0.5], df[[i]], lower.tail=F)
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
		subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
		# cat(sprintf("Number of sub-threshold F p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvaluesNumbers) >= 1 ) {
			thresholds[j,i] <- qf(max(subTholdPvaluesNumbers), df[[i]][1], df[[i]][2], lower.tail=FALSE)
		} else { thresholds[j,i] <- NA }
      }
      else if (statType[i] %in% c("t", "tlmer")) {
		subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
		subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
		#cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvaluesNumbers) >= 1 ) {
			thresholds[j,i] <-qt(max(subTholdPvaluesNumbers)/2, df[[i]], lower.tail=FALSE)
		} else { thresholds[j,i] <- NA }
      }
      else if (statType[i] == "u") {
		subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
		subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
		#cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvaluesNumbers) >= 1 ) {
			thresholds[j,i] <-qwilcox(max(subTholdPvaluesNumbers),m,n,lower.tail = TRUE)
		} else { thresholds[j,i] <- NA }
      }

      else if (statType[i] == "chisq") {
 		subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
		#cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
		if ( length(subTholdPvalues) >= 1 ) {
			thresholds[j,i] <-qchisq(max(subTholdPvalues), df[[i]], lower.tail=FALSE)
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
###########################################################################################


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

mincTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="t-test", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("t") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  rttest <- t.test(vox~grouping,data=gf,paired=FALSE)
  attr(result, "df") <- rttest$parameter
  colnames(result) <- c("T-statistic")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)
}

mincPairedTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="paired-t-test", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("t") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  rttest <- t.test(vox~grouping,data=gf,paired=TRUE)
  attr(result, "df") <- rttest$parameter
  colnames(result) <- c("T-statistic")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)
}

mincCorrelation <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="correlation", maskval=maskval)
  return(result)
}

mincWilcoxon <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="wilcoxon", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("u") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  a = levels(gf$grouping)
  attr(result, "m") <- nrow(gf[gf$grouping == a[1],])
  attr(result, "n") <- nrow(gf[gf$grouping == a[2],])
  colnames(result) <- c("Mann-Whitney")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)



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
                  matrix(),
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
    if(!grepl("t-test",method) && !grepl("correlation",method) && !grepl("wilcoxon",method)) {
      colnames(result) <- levels(grouping)
    }
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
                       mask=NULL, workers=4, tinyMask=FALSE, method="snowfall",global="",packages="") {
  
  REDUCE = TRUE; # For now this option is not exposed

  # if no mask exists use the entire volume
  if (is.null(mask)) {
    maskV = mincGetVolume(filenames[1])
    nVoxels = length(maskV)
    maskV[maskV >= min(maskV)] <- as.integer(cut(seq_len(nVoxels), workers)) 
  }
  else {
    maskV <- mincGetVolume(mask)
    # optionally make the mask a fraction of the original size - for testing
    if (tinyMask!=FALSE) {
      maskV[maskV>1.5] <- 0
    }
    nVoxels <- sum(maskV>0.5)
    maskV[maskV>0.5] <- as.integer(cut(seq_len(nVoxels), workers)) 
  }
  
  maskFilename <- paste("pmincApplyTmpMask-", Sys.getpid(), ".mnc", sep="")
  mincWriteVolume(maskV, maskFilename, clobber=TRUE)
  
  pout <- list()
  
  if (method == "local") {
    stop("Lovely code ... that generates inconsistent results because something somewhere is not thread safe ...")
    
    library(multicore)
    library(doMC)
    library(foreach)
    registerDoMC(workers)
    
    # run the job spread across each core
    pout <- foreach(i=1:workers) %dopar% { mincApply(filenames, function.string,
                                                   mask=maskFilename, maskval=i) }
    #cat("length: ", length(pout), "\n")
  }
  else if (method == "sge") {
    library(Rsge)

	options(sge.use.cluster = TRUE)
	options(sge.block.size = 100)
	options(sge.user.options= "-S /bin/bash")
	options(sge.ret.ext= "sge.ret")
	options(sge.use.qacct= FALSE)
	options(sge.trace = TRUE) 
	options(sge.save.global = FALSE)
	options(sge.qsub.options = "-cwd")
	options(sge.qsub.blocking = "-sync y -t 1-")
	options(sge.monitor.script = "MonitorJob.sh")
	options(sge.script="RunSgeJob")
	options(sge.file.prefix="Rsge_data")
	options(sge.debug=TRUE)
	options(sge.remove.files=FALSE)
	options(sge.qacct= "qacct")
	options(sge.qstat= "qstat")
	options(sge.qsub= "qsub")

    # Need to use double quotes, because both sge.submit and mincApply try to evalute the functin
    function.string = enquote(function.string)
    
    l1 <- list(length=workers)
    
    if(packages == "")
	packageList = c("RMINC")
    else
        packageList=c(packages,"RMINC")

    if(global == "") 
	globallist = c(sub("\\(([A-Z]|[a-z])\\)","",function.string))
    else
        globallist = c(global,sub("\\(([A-Z]|[a-z])\\)","",function.string))

    # Submit one job to the queue for each segmented brain region
    for(i in 1:workers) {
      l1[[i]]<- sge.submit(mincApply,filenames,function.string, mask=maskFilename,reduce=TRUE,
                           maskval=i, packages=packageList,global.savelist= globallist)
      
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
    library(snowfall)
    sfInit(parallel=TRUE, cpus=workers)
    if(packages == "")
	packageList = c("RMINC")
    else
        packageList=c(packages,"RMINC")

    for (nPackage in 1:length(packageList)) {
    	sfLibrary(packageList[nPackage],character.only=TRUE) }


    sfExport(list = global) 
 
    wrapper <- function(i) {
      cat( "Current index: ", i, "\n" ) 
      return(mincApply(filenames, function.string, mask=maskFilename,
                       maskval=i, reduce=REDUCE))
    }
    # use all workers in the current cluster if # of workers not specified
    if (is.null(workers)) {
      workers <- length(sfSocketHosts())
    }
    
   

    #sink("/dev/null");
    pout <- sfLapply(1:workers, wrapper)
    
    sfStop();
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
  
  for(i in 1:workers) {
    if (length(test)>1) {
      if(REDUCE == TRUE)	
      	output[maskV == i,] <- pout[[i]]
      else
      	output[maskV == i,] <- pout[[i]][maskV == i, ] 
    }
    else {
      output[maskV==i] <- pout[[i]]
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
              matrix(),
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
		   matrix(),
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
vertexAnova <- function(formula, data, subset=NULL) {
  # Create Model
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mmatrix <- model.matrix(formula, mf)

  # Load Vertex Data from Files
  filenames <- as.character(mf[,1])
  data.matrix <- vertexTable(filenames)
  result <- .Call("vertex_anova_loop", data.matrix, mmatrix,attr(mmatrix, "assign"), PACKAGE="RMINC");

  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- rep("F", ncol(result))

  # Use the first voxel in order to get the dimension names
  v.firstVoxel <- data.matrix[1,]
  
  # Get the degrees of freedom
  l <- lm.fit(mmatrix, v.firstVoxel)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
        dflist[[i]] <- c(dfs[[i + 1]], dfr)
    }
  attr(result, "df") <- dflist
  
  colnames(result) <- attr(terms(formula), "term.labels")
 
  class(result) <- c("vertexMultiDim", "matrix")
  return(result)
}
###########################################################################################
#' Calculates statistics and coefficients for linear model of specified vertex files
#' @param formula a model formula. The RHS of the formula may contain one term with filenames. If
#' so only the + operator may be used, and only two terms may appear on the RHS
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

  # Build model.frame
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  

  if(length(grep("\\$",formula[[3]])) > 0) {
    stop("$ Not Permitted in Formula")  
  }

  # the following returns:
  #
  # list(data.matrix.left = data.matrix.left, 
  #      data.matrix.right = data.matrix.right,
  #      rows = rows,
  #      matrixFound = matrixFound,
  #      mmatrix = mmatrix)
  parseLmOutput <- parseLmFormula(formula,data,mf)  

  if(parseLmOutput$matrixFound) {
    parseLmOutput$data.matrix.left  <- vertexTable(parseLmOutput$data.matrix.left)
    parseLmOutput$data.matrix.right <- vertexTable(parseLmOutput$data.matrix.right)
  }
  else {
    filenames <- as.character(mf[,1])
    parseLmOutput$mmatrix <- model.matrix(formula, mf)	
    parseLmOutput$data.matrix.left <- vertexTable(filenames)
    parseLmOutput$rows = colnames(parseLmOutput$mmatrix)
  } 
  
  result <- .Call("vertex_lm_loop",
                  parseLmOutput$data.matrix.left,
                  parseLmOutput$data.matrix.right,
                  parseLmOutput$mmatrix,
                  PACKAGE="RMINC") 

  attr(result, "likeVolume") <- as.character(mf[,1])[1]
  attr(result, "model") <- as.matrix(parseLmOutput$mmatrix)
  attr(result, "filenames") <- as.character(mf[,1])
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
 
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2+1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  betaNames = paste('beta-', parseLmOutput$rows, sep='')
  tnames = paste('tvalue-', parseLmOutput$rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("vertexMultiDim", "matrix")

  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

vertexMean <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(rowMeans(vertexData))

} 


vertexSum <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(rowSums(vertexData))

} 

vertexVar <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(apply(vertexData,1,var))

} 

vertexSD<- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(apply(vertexData,1,sd))

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

parseLmFormula <- function(formula,data,mf) 
{
  matrixName = ''
  mmatrix = matrix()
  data.matrix.right = matrix()
  data.matrix.left = matrix()
  rows = NULL
  matrixFound = FALSE
  
  if(length(formula[[3]]) == 1) {
    # Only 1 Term on the RHS
    # If the formula has the form LHS ~ RHS (for instance jacobians ~ Genotype)
    # it will get split up as follows:
    # formula[[1]] -> ~
    # formula[[2]] -> LHS (jacobians)
    # formula[[3]] -> RHS (Genotype)
    #
    # It is possible to have something like: jacobians ~ Genotype * Age
    # where the RHS (formula[[3]]) would be "Genotype * Age". Here, we 
    # deal with cases where there is only a single term on the right hand side
    if(is.null(data)) {
      # This was a bug discovered in the end of July 2014. The general way
      # we call mincLm is by using a dataframe that contains all terms. So 
      # it could be: mincLm(jacobians ~ Genotype, gf). But it's quite
      # possible not to have this dataframe around. The command in the 
      # else block was the standard way to retrieve the right hand side
      # of the formula, this if block fixes the case where there is no dataframe
      # provided.
      rCommand = paste("term <-",formula[[3]],sep="")
    }
    else {
      rCommand = paste("term <- data$",formula[[3]],sep="")
    }
    
    eval(parse(text=rCommand))
    fileinfo = file.info(as.character(term[1]))
    if (!is.na(fileinfo$size)) {
      # Save term name for later
      rows = c('Intercept',formula[[3]])
      matrixName = formula[[3]]
      matrixFound = TRUE
      data.matrix.left <- as.character(mf[,1])
      data.matrix.right <- as.character(mf[,2])
    }
  }
  # Multiple Terms on RHS
  else {
    for (nTerm in 2:length(formula[[3]])){
      # Skip if it is an interaction term
      if(length(formula[[3]][[nTerm]] > 1))
        next
      rCommand = paste("term <- data$",formula[[3]][[nTerm]],sep="")
      # Skip if it is a formula symbol (i.e. *)
      if(!as.character(formula[[3]][[nTerm]]) %in% names(data)) {
        next
      }
      eval(parse(text=rCommand))	
      fileinfo = file.info(as.character(term[1]))
      if (!is.na(fileinfo$size)) {
        if(length(grep('\\+',formula[[3]][[1]])) == 0) {
          stop("Only + sign allowed when using filenames")
        }
        if(length(formula[[3]]) > 3) {
          stop("Only 2 terms allowed when using filenames on the RHS")	
        }
        matrixName = formula[[3]][[nTerm]]
        matrixFound = TRUE
        data.matrix.left <- as.character(mf[,1])
        data.matrix.right <- as.character(mf[,nTerm])
      }
      else {
        tmpFormula = formula
        rCommand = paste("formula <-",formula[[2]],"~",formula[[3]][[nTerm]],sep="")
        eval(parse(text=rCommand))	
        mmatrix <- model.matrix(formula, mf)	
        formula = tmpFormula	
      }
    }
    rows = colnames(mmatrix)
    rows = append(rows,matrixName)
  }
  return(list(data.matrix.left = data.matrix.left, data.matrix.right = data.matrix.right,rows = rows,matrixFound = matrixFound,mmatrix = mmatrix))
}

### lmer functions stuff starts here
#' mincified version of lmer from lme4
#'
#' mincLmer should be used the same way as a straight lmer call, except
#' that the left hand side of the equation contains minc filenames rather than
#' an actual response variable.
#'
#' @param formula the lmer formula, filenames go on left hand side
#' @param data the data frame, all items in formula should be in here
#' @param mask the mask within which lmer is solved
#' @param parallel how many processors to run on (default=single processor).
#' Specified as a two element vector, with the first element corresponding to
#' the type of parallelization (sge or snowfall), and the second to the number
#' of processors to use. 
#' @param REML whether to use use Restricted Maximum Likelihood or Maximum Likelihood
#' @param control lmer control function
#' @param start lmer start function
#' @param verbose lmer verbosity control
#'
#' @return a matrix where rows correspond to number of voxels in the file and columns to
#' the number of terms in the formula, with both the beta coefficient and the t-statistic
#' being returned. In addition, an extra column keeps the log likelihood, and another
#' whether the mixed effects fitting converged or not.
#'
#' @seealso \code{\link{lmer}} for description of lmer and lmer formulas; \code{\link{mincLm}}
#'
#' @examples
#' \dontrun{
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc")
#' mincWriteVolume(vs, "age-term.mnc", "tvalue-age")
#' # run in parallel with multiple processors on the local machine
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc", parallel=c("snowfall", 4))
#' # run in parallel with multiple processors over the sge batch queueing system
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc", parallel=c("sge", 4))
#' }
mincLmer <- function(formula, data, mask=NULL, parallel=NULL,
                     REML=TRUE, control=lmerControl(), start=NULL, verbose=0L) {
  # the outside part of the loop - setting up various matrices, etc., whatever that is
  # constant for all voxels goes here

  # code ripped straight from lme4::lmer
  mc <- mcout <- match.call()
  mc$control <- lmerControl()
  mc[[1]] <- quote(lme4::lFormula)

  # remove mask and parallel, since lmer does not know about them and keeping them
  # generates obscure warning messages
  index <- match("mask", names(mc))
  if (!is.na(index)) {
    mc <- mc[-index]
  }
  index <- match("parallel", names(mc))
  if (!is.na(index)) {
    mc <- mc[-index]
  }
  lmod <- eval(mc, parent.frame(1L))

  # code ripped from lme4:::mkLmerDevFun
  rho <- new.env(parent = parent.env(environment()))
  rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta", 
                                                  "Lambdat", "Lind")],
                                    n = nrow(lmod$X), list(X = lmod$X)))
  REMLpass <- if (REML) 
    ncol(lmod$X)
  else 0L


  mincLmerList <<- list(lmod, mcout, control, start, verbose, rho, REMLpass)

  # for some reason there is a namespace issue if I call diag directly, but only if inside
  # a function that is part of RMINC (i.e. if I source the code it works fine). So here's a
  # workaround to get the method first, give it a new name, and assign to global namespace.
  tmpDiag <<- getMethod("diag", "dsyMatrix")
  #fmincLmerOptimizeAndExtract <<- mincLmerOptimizeAndExtract

  if (!is.null(parallel)) {
    # a vector with two elements: the methods followed by the # of workers
    if (parallel[1] %in% c("local", "snowfall", "sge")) {
      if (parallel[1] == "local") {
        parallel[1] <- "snowfall" #local and snowfall are synonymous
      }
      out <- pMincApply(lmod$fr[,1],
                        quote(mincLmerOptimizeAndExtract(x)),
                        mask=mask,
                        method=parallel[1],
                        worker=as.numeric(parallel[2]),
                        global=c("mincLmerList", "tmpDiag"),# "mincLmerOptimize", "mincLmerExtractVariables"),
                        packages=c("lme4", "RMINC"))
    }
    else {
      stop("Error: unknown parallelization method")
    }
  }
  else {
    out <- mincApply(lmod$fr[,1], # assumes that the formula was e.g. filenames ~ effects
                     quote(mincLmerOptimizeAndExtract(x)),
                     mask=mask)
  }

  # set Inf to 0 (Inf's are generated when vcov can't compute)
  out[is.infinite(out)] <- 0
  
  termnames <- colnames(lmod$X)
  betaNames <- paste("beta-", termnames, sep="")
  tnames <- paste("tvalue-", termnames, sep="")
  colnames(out) <- c(betaNames, tnames, "logLik", "converged")

  # generate some random numbers for a single fit in order to extract some extra info
  mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])))
  
  attr(out, "stat-type") <- c(rep("beta", length(betaNames)), rep("tlmer", length(tnames)),
                              "logLik", "converged")
  # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
  attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
  attr(out, "REML") <- REML
  attr(out, "mask") <- mask
  attr(out, "mincLmerList") <- mincLmerList
  class(out) <- c("mincLmer", "mincMultiDim", "matrix")

  return(out)
}

#' estimate the degrees of freedom for parameters in a mincLmer model
#'
#' There is much uncertainty in how to compute p-values for mixed-effects
#' statistics, related to the correct calculation of the degrees of freedom
#' of the model (see here \link{http://glmm.wikidot.com/faq#df}). mincLmer by
#' default does not return the degrees of freedom as part of its model, instead
#' requiring an explicit call to a separate function (such as this one).
#' The implementation here is the Satterthwaite approximation. This approximation
#' is computed from the data, to avoid the significant run-time requirement of computing
#' it separate for every voxel, here it is only computed on a small number of voxels
#' within the mask and the median DF returned for every variable.
#' 
#' @param model the output of mincLmer
#'
#' @return the same mincLmer model, now with degrees of freedom set
#'
#' @seealso \code{\link{mincLmer}} for mixed effects modelling, \code{\link{mincFDR}}
#' for multiple comparisons corrections.
#'
#' @examples
#' \dontrun{
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc")
#' vs <- mincLmerEstimateDF(vs)
#' qvals <- mincFDR(vs, mask=attr(vs, "mask"))
#' qvals
#' }
mincLmerEstimateDF <- function(model) {
  # set the DF based on the Satterthwaite approximation
  # load lmerTest library if not loaded; lmerTest takes over some lmer functions, so unload if
  # it wasn't loaded in the first place
  lmerTestLoaded <- "package:lmerTest" %in% search()
  if (lmerTestLoaded == FALSE) {
    library(lmerTest)
  }

  # put the lmod variable back in the global environment
  #lmod <<- attr(model, "mincLmerList")[[1]]
  mincLmerList <<- attr(model, "mincLmerList")
  mask <- attr(model, "mask")
  
  # estimated DF depends on the input data. Rather than estimate separately at every voxel,
  # instead select a small number of voxels and estimate DF for those voxels, then keep the
  # min
  nvoxels <- 50
  rvoxels <- mincSelectRandomVoxels(mask, nvoxels)
  dfs <- matrix(nrow=nvoxels, ncol=sum(attr(model, "stat-type") %in% "tlmer"))
  for (i in 1:nvoxels) {
    voxelData <- mincGetVoxel(mincLmerList[[1]]$fr[,1], rvoxels[i,])
    mmod <- mincLmerOptimize(voxelData)
    # code directly from lmerTest library
    rho <- lmerTest:::rhoInit(mmod)
    hessian <- lmerTest:::myhess
    Dev <- lmerTest:::Dev
    h  <-  hessian(function(x) Dev(rho,x), rho$param$vec.matr)
    rho$A <- 2*solve(h)
    dfs[i,] <- lmerTest:::calculateTtest(rho, diag(rep(1, length(rho$fixEffs))),
                                         length(rho$fixEffs), "simple")[,1]
    
  }
  df <- apply(dfs, 2, median)
  cat("Mean df: ", apply(dfs, 2, mean), "\n")
  cat("Median df: ", apply(dfs, 2, median), "\n")
  cat("Min df: ", apply(dfs, 2, min), "\n")
  cat("Max df: ", apply(dfs, 2, max), "\n")
  cat("Sd df: ", apply(dfs, 2, sd), "\n")
  
  attr(model, "df") <- df
  if (lmerTestLoaded == FALSE) {
    detach("package:lmerTest")
  }
  return(model)
}
# the actual optimization of the mixed effects models; everything that has to be recomputed
# for every voxel goes here. Works on x (each voxel is assigned x during the loop), and
# assumes that all the other info is in a variable called mincLmerList in the global
# environment. This last part is a hack to get around the lack of multiple function arguments
# for mincApply and friends.
mincLmerOptimize <- function(x) {
  # code ripped straight from lme4::lmer
  # assignments from global variable set in mincLmer
  lmod <- mincLmerList[[1]]
  mcout <- mincLmerList[[2]]
  control <- mincLmerList[[3]]
  start <- mincLmerList[[4]]
  verbose <- mincLmerList[[5]]
  rho <- mincLmerList[[6]]
  REMLpass <- mincLmerList[[7]]
  
  # assign the vector of voxel values
  lmod$fr[,1] <- x

  mmod <- mincLmerOptimizeCore(rho, lmod, REMLpass, verbose, control, mcout, start, reinit=FALSE)
  # the rho$pp merModP object needs to be reinitialized at (to me) mysterious times. That is
  # quite time consuming, however, so only do so if something did not converge.
  if (length(attr(mmod,"optinfo")$conv$lme4$code) != 0) {
    #cat("Restarting ... \n")
    mmod <- mincLmerOptimizeCore(rho, lmod, REMLpass, verbose, control, mcout, start, reinit=TRUE)
  }
  
  return(mmod)
}

# the core code that does the optimization. Only reason it is not part of mincLmerOptimize
# is to allow for the reinitializtion in case of convergence error
mincLmerOptimizeCore <- function(rho, lmod, REMLpass, verbose, control, mcout, start, reinit=F) {
    # finish building the dev function by adding the response term
  # code from lme4:::mkLmerDevFun
  # this code is quite slow; I'm baffled about when it's needed, as most of the
  # time it can stay outside the loop, but occasionally gives weird errors
  # if inside. So wrapped inside that reinit bit:
  if (reinit) {
    lmod$reTrms <- mkReTrms(findbars(lme4:::RHSForm(mcout[[2]])), lmod$fr)
    rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta", 
                                                    "Lambdat", "Lind")],
                                      n = nrow(lmod$X), list(X = lmod$X)))
  }

  rho$resp <- mkRespMod(lmod$fr, REML = REMLpass)
  devfun <- lme4:::mkdevfun(rho, 0L, verbose, control)
  theta <- lme4:::getStart(lmod$start, lmod$reTrms$lower, rho$pp)
  if (length(rho$resp$y) > 0) 
    devfun(rho$pp$theta)
  rho$lower <- lmod$reTrms$lower

  # kept the old full mkLmerDevFun call around here in case the divided call
  # ends up with unexpected side effects down the road.
  #devfun <- do.call(mkLmerDevfun, c(lmod,
  #                                  list(start = start, 
  #                                       verbose = verbose,
  #                                       control = control)))
  #devfun <- mkLmerDevfun(lmod$fr, lmod$X, lmod$reTrms, lmod$REML, start, verbose, control)

  # the optimization of the function - straight from lme4:::lmer
  opt <- optimizeLmer(devfun, optimizer = control$optimizer, 
                      restart_edge = control$restart_edge,
                      boundary.tol = control$boundary.tol, 
                      control = control$optCtrl, verbose = verbose,
                      start = start, 
                      calc.derivs = control$calc.derivs,
                      use.last.params = control$use.last.params)
  cc <- lme4:::checkConv(attr(opt, "derivs"), opt$par,
                         ctrl = control$checkConv, 
                         lbound = environment(devfun)$lower)
  mmod <- mkMerMod(environment(devfun), opt, lmod$reTrms,
                   fr = lmod$fr, mcout, lme4conv = cc)
  return(mmod)
}

# takes a merMod object, gets beta, t, and logLikelihood values, and
# returns them as a vector
mincLmerExtractVariables <- function(mmod) {
  se <- tryCatch({ # vcov sometimes complains that matris is not positive definite
    sqrt(tmpDiag(vcov(mmod, T)))
  }, warning=function(w) {
    return(0)
  }, error=function(e) {
    return(0)
  })
  fe <- mmod@beta
  t <- fe / se # returns Inf if se=0 (see error handling above); mincLmer removes Inf
  ll <- logLik(mmod)
  # the convergence return value (I think; need to verify) - should be 0 if it converged
  converged <- length(attr(mmod,"optinfo")$conv$lme4$code) == 0
  return(c(fe,t, ll, converged))
}

mincLmerOptimizeAndExtract <- function(x) {
  mmod <- mincLmerOptimize(x)
  return(mincLmerExtractVariables(mmod))
}

#' run log likelihood ratio tests for different mincLmer objects
#'
#' Computes the log likelihood ratio of 2 or more voxel-wise lmer calls, testing the hypothesis that
#' the more complex model better fits the data. Note that it requires the mixed effects to have been
#' fitted with maximum likelihood, and not restricted maximum likelihood; in other words, if you want
#' to use these log likelihood tests, make sure to specify REML=FALSE in mincLmer.
#'
#' @return the voxel wise log likelihood test. Will have a number of columns corresponding to the
#' number of inputs -1. Note that it resorts the inputs from lowest to highest degrees of freedom
#'
#' @seealso \code{\link{lmer}} and \code{\link{mincLmer}} for description of lmer and mincLmer.
#' \code{\link{mincFDR}} for using the False Discovery Rate to correct for multiple comparisons,
#' and \code{\link{mincWriteVolume}} for outputting the values to MINC files.
#'
#' @examples
#' \dontrun{
#' m1 <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc", REML=F)
#' m2 <- mincLmer(filenames ~ age + I(age^2) + sex + (age|id), data=gf, mask="mask.mnc", REML=F)
#' m3 <- mincLmer(filenames ~ age + I(age^2) + I(age^3) + sex + (age|id), data=gf, mask="mask.mnc", REML=F)
#' llr <- mincLogLikRatio(m1, m2, m3)
#' mincFDR(llr)
#' mincWriteVolume(llr, "m2vsm3.mnc", "m3")
#' }
mincLogLikRatio <- function(...) {
  dots <- list(...)

  # get the names of the actual objects passed it; used for naming output columns
  dotslist <- substitute(list(...))[-1];
  inputnames <- sapply(dotslist, deparse) 

  # test for REML vs ML, exit if REML. 
  for (i in 1:length(dots)) {
    REML <- attr(dots[[i]], "REML")
    if (is.null(REML)) {
      stop("all arguments must be the outputs of mincLmer")
    }
    else if (REML == TRUE) {
      stop("Log likelihood ratio tests only work reliably if fitted with maximum likelihood, but not if fitted with restricted maximum likelihood. Rerun your model with REML=FALSE")
    }
  }

  # sort the degrees of freedom of each of the models from lowest to highest
  df <- vector(length=length(dots))
  for (i in 1:length(dots)) {
    df[i] <- attr(dots[[i]], "logLikDF")
  }
  dfs <- sort(df, index.return=T)
  # create the output matrix - number of columns equal to number of objects passed in minus 1
  logLikMatrix <- matrix(nrow=nrow(dots[[1]]), ncol=length(dots))
  for (i in 1:length(dots)) {
    logLikMatrix[,i] <- dots[[dfs$ix[i]]][,"logLik"]
  }
                             
  # compute the log likelihood ratio test
  flogLikRatio <- function(x) { 2 * pmax(0, diff(x)) }
  # apply at every voxel
  # erm - apply is slow as a dog here ... replace with manual rolled function
  #out <- t(apply(logLikMatrix, 1, flogLikRatio))
  out <- matrix(nrow=nrow(dots[[1]]), ncol=length(dots)-1)
  outnames <- vector(length=length(dots)-1)
  for (i in 2:length(dots)) {
    out[, i-1] <- 2 * abs(logLikMatrix[, i] - logLikMatrix[, i-1])
    outnames[i-1]  <- inputnames[dfs$ix[i]]
  }

  # set attributes and class types
  attr(out, "likeVolume") <- attr(dots[[1]], "likeVolume")
  attr(out, "stat-type") <- rep("chisq", ncol(out))
  attr(out, "df") <- diff(dfs$x)
  attr(out, "mask") <- attr(dots[[1]], "mask")

  # keep the mincLmerList for every object
  mincLmerLists <- list()
  for (i in 1:length(dots)) {
    mincLmerLists[[i]] <- attr(dots[[dfs$ix[i]]], "mincLmerList")
  }
  attr(out, "mincLmerLists") <- mincLmerLists
  if (length(dots) == 2) {
    class(out) <- c("mincLogLikRatio", "mincSingleDim", "numeric")
  }
  else {
    class(out) <- c("mincLogLikRatio", "mincMultiDim", "matrix")
  }
  colnames(out) <- outnames
  return(out)
}

#' computes parametric bootstrap on mincLogLikRatio output
#'
#' The Log Likelihood Ratio tests closely approximates a Chi-squared
#' distribution when the number of groups (i.e. individual subjects in a
#' longitudinal study) is large (>50), but can be anticonservative when
#' small. A parametric bootstrap test, in which data is randomly simulated from
#' the null model and then fit with both models, can give the correct p-value.
#' Here we compute the parametric boostrap on a small number of randomly chosen
#' voxels to get a sense of biased the estimated p-values from the log likelihood
#' ratio test really were.
#'
#' @param logLikOutput the output from mincLogLikRatio
#' @param selection the algorithm for randomly chosing voxels. Only "random" works for now.
#' @param nsims the number of simulations to run per voxel
#' @param nvoxels the number of voxels to run the parametric bootstrap on
#'
#' @return a matrix containing the chi-square p-values and the bootstrapped p-values
mincLogLikRatioParametricBootstrap <- function(logLikOutput, selection="random",
                                               nsims=500, nvoxels=50) {
  mincLmerLists <- attr(logLikOutput, "mincLmerLists")
  mask <- attr(logLikOutput, "mask")
  if (length(mincLmerLists) != 2) {
    stop("Error: parametric bootstrap only implemented for the two model comparison case")
  }
  if (selection == "random") {
    out <- matrix(nrow=nvoxels, ncol=2)
    simLogLik <- matrix(nrow=nsims, ncol=2)
    simLogLikRatio <- numeric(nvoxels)
    voxelMatrix <- matrix(nrow=nsims, ncol=nrow(attr(logLikOutput, "mincLmerLists")[[1]][[1]]$X))
    voxels <- mincSelectRandomVoxels(mask, nvoxels, convert=F)
    #cat(voxels)
    for (i in 1:nvoxels) {
      # refit the null model first since we'll need the merMod object for simulations
      mincLmerList <<- mincLmerLists[[1]]
      voxel <- mincGetVoxel(mincLmerList[[1]]$fr[,1], mincVectorToVoxelCoordinates(mask, voxels[i]))
      mmod <- mincLmerOptimize(voxel)
      for (j in 1:nsims) {
        # create the simulated data from the null model
        voxelMatrix[j,] <- unlist(simulate(mmod))
        # compute the log likelihood for the null model
        simLogLik[j,1] <- logLik(mincLmerOptimize(voxelMatrix[j,]))
      }
      # do it all again for the alternate model (happens in separate loop since
      # mincLmerOptimize relies on the global variable mincLmerList)
      mincLmerList <<- mincLmerLists[[2]]
      for (j in 1:nsims) {
        simLogLik[j,2] <- logLik(mincLmerOptimize(voxelMatrix[j,]))
      }
      # compute the normally estimated chisq p value
      out[i,1] <- pchisq(logLikOutput[voxels[i]], attr(logLikOutput, "df"), lower=F)
      # compute the parametric log likelihood ratio
      simLogLikRatio <- 2 * abs(simLogLik[,1] - simLogLik[,2])
      # compute the parametric bootstrap p value
      out[i,2] <- mean( simLogLikRatio >= logLikOutput[voxels[i]] )
      colnames(out) <- c("chisq", "parametricBootstrap")
    }
  }
  else {
    stop("Error: unknown voxel selection mechanism")
  }
  return(out)
}

  
### end of lmer bits of code

#' converts a vector index to the voxel indices in MINC
#'
#' RMINC stores volume data as 1d arrays. This function gives the
#' corresponding voxel coordinates (in the dimension order of the volume)
#' for an index into the 1d array.
#'
#' @param volumeFileName the filename of the MINC volume
#' @param vectorCoord the integer array index to convert
#'
#' @return a vector of length 3 containing the MINC indices in volume dimension order
#'
#' @examples
#' \dontrun{
#' index <- mincVectorToVoxelCoordinates("filename.mnc", 345322)
#' voxel <- mincGetVoxel(gf$filenames, index)
#' }
mincVectorToVoxelCoordinates <- function(volumeFileName, vectorCoord) {
  sizes <- minc.dimensions.sizes(volumeFileName)
  # the fun off by one bit
  vectorCoord <- vectorCoord-1
  i1 <- vectorCoord %/% (sizes[2]*sizes[3])
  i1r <- vectorCoord %% (sizes[2]*sizes[3])
  i2 <- i1r %/% sizes[2]
  i3 <- i1r %% sizes[2]
  return(c(i1, i2, i3))
}

#' selects a few random indices from a volume
#'
#' Given a filename, select a few random indices using the uniform distribution
#' from voxels that have a value of 1 (i.e. from a mask volume)
#'
#' @param volumeFileName the filename for a MINC volume
#' @param nvoxel the number of voxels to select
#' @param convert whether to convert to MINC voxel space (default) or keep in index space  
mincSelectRandomVoxels <- function(volumeFileName, nvoxels=50, convert=TRUE) {
  #load volume - should be a binary mask
  mvol <- mincGetVolume(volumeFileName) 
  # get the indices of voxels inside mask
  vinmask <- which(mvol %in% 1)
  # keep a random set of voxels
  indicesToKeep <-  vinmask[ floor(runif(nvoxels, min=1, max=length(vinmask))) ]
  if (convert == TRUE) {
    out <- matrix(nrow=nvoxels, ncol=3)
    for (i in 1:nvoxels) {
      out[i,] <- mincVectorToVoxelCoordinates(volumeFileName, indicesToKeep[i])
    }
    return(out)
  }
  else {
    return(indicesToKeep)
  }
}

# Run Testbed
runRMINCTestbed <- function() {

  # Make sure environment is clear
  #rm(list=ls())

  system('mkdir /tmp/rminctestdata')


  # Download Tarball from Wiki
  system("wget -O /tmp/rminctestdata/rminctestdata.tar.gz --no-check-certificate https://wiki.phenogenomics.ca/download/attachments/1654/rminctestdata.tar.gz")

  # Untar
  system('tar -xf /tmp/rminctestdata/rminctestdata.tar.gz -C /tmp/')
  library(testthat)

  # Run Tests
  rmincPath = find.package("RMINC")
  cat("\n\nRunning tests in: ", paste(rmincPath,"/","tests/",sep=""), "\n\n\n")
  test_dir(paste(rmincPath,"/","tests/",sep=""))
  
  cat("\n*********************************************\n")
  cat("The RMINC test bed finished running all tests\n")
  cat("*********************************************\n\n\n")
  # Remove temp data, and downloaded files
  cat("Removing temporary directory /tmp/rminctestdata\n")
  system('rm -fr /tmp/rminctestdata')
  
}


