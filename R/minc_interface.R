#' Set Masked Value
#' 
#' Use this function to configure how masked values should be treated
#' @param val The new value for RMINC_MASKED_VALUE, defaults to resetting
#' the masked value to NA
#' @export
setRMINCMaskedValue <- 
  function(val = NA){
    options(RMINC_MASKED_VALUE = structure(val, class = "RMINC_MASKED_VALUE"))
    return(invisible(NULL))
  }

#' @method print RMINC_MASKED_VALUE
print.RMINC_MASKED_VALUE <-
  function(x, ...)
    print(as.symbol("masked"))
    
mincGetTagFile <- function(filename) {
  tags <- scan(filename, what="character")
  # tag points begin after Points = line
  beginIndex <- grep("Points", tags) + 2
  endIndex <- length(tags)-1 # get rid of training ;
  return(matrix(as.numeric(tags[beginIndex:endIndex]), ncol=7, byrow=T))
}

mincConvertTagToMincArrayCoordinates <- function(tags, filename) {
  tags <- tags[,c(1:3,7)] # get rid of repeated coordinates
  out <- matrix(ncol=ncol(tags), nrow=nrow(tags))
  for (i in 1:nrow(tags)) {
    out[i,3:1] <- mincConvertWorldToVoxel(filename, tags[i,1], tags[i,2], tags[i,3]) + 1
  }
  out[,4] <- tags[,4]
  return(out)
}

#' Retrieve Voxel Values
#' 
#' Return the intensity of a given voxel in a set of minc files
#' 
#' @param filenames paths to the minc files
#' @param v1 Either a 3-element vector of voxel coordinates or the first
#' @param v2 the second voxel coordinate if not NULL
#' @param v3 the third voxel coordinate if not NULL
#' @return Returns a \code{mincVoxel} object containing a vector
#' of intensities and attributes specify the voxel and world coordinates
#' of the values.
#' @export
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

#' Minc File Check
#' 
#' test to see whether files exist and are readable
#' @param filenames paths to the files to check
#' @details Throws an error if a file is not found or not readable
#' @return Returns NULL invisibly
mincFileCheck <- function(filenames) {
  for(i in 1:length(filenames) ) {
    if(file.access(as.character(filenames[i]), 4) == -1 ){
        stop("The following file could not be read (full filename is between the dashes): ---", filenames[i], "---")
    }
  }
  
  return(invisible(NULL))
}


#' World Vector
#' 
#' Given a set of world coordinates return the value at those coordinates
#' from a set of minc files as a numeric vector
#' 
#' @param filenames A character vector of one or more filenames
#' @param v1 First world coordinate
#' @param v2 Second world coordinate
#' @param v3 Third world coordinate
#' @return a vector of values 
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

#' Voxel Vector
#' 
#' Given a voxel coordinate, return the value at those coordinates
#' from a set of minc files as a numeric vector
#' 
#' @param filenames Character vector with paths to minc files
#' @param v1 First voxel coordinate
#' @param v2 Second voxel coordinate
#' @param v3 Third voxel coordinate
#' @param v.length The number of values to return
mincGetVector <- function(filenames, v1, v2, v3, v.length = NULL) {
  num.files <- length(filenames)
  if(is.null(v.length)) v.length <- num.files
  
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

#' Voxel to World
#' 
#' Convert a discrete set of voxel coordinates into the volumes
#' continous world coordinate space
#' @param filename Character vector with paths to minc files
#' @param v1 First voxel coordinate
#' @param v2 Second voxel coordinate
#' @param v3 Third voxel coordinate
#' @return a 3-component numeric vector of world coordinates
#' @export
mincConvertVoxelToWorld <- function(filename, v1, v2, v3) {
  output <- .C("convert_voxel_to_world",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3), PACKAGE="RMINC")$o
  return(output)
}

#' World to Voxel
#' 
#' Convert coordinates in the volumes world space (a continuous coordinate space)
#' to a discrete set of voxel coordinates
#' @param filename A path to a minc file
#' @param v1 First world coordinate
#' @param v2 Second world coordinate
#' @param v3 Third world coordinate
#' @return a 3-component numeric vector of voxel coordinates
#' @export
mincConvertWorldToVoxel <- function(filename, v1, v2, v3) {
  output <- .C("convert_world_to_voxel",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3), PACKAGE="RMINC")$o
  return(round(output))
}

#' @title Read a MINC file
#' @description Load a 3-dimensional MINC2 volume and returns it as a 1D array
#' @name mincGetVolume
#' @title Return a volume as a 1D array
#' @param filename A string corresponding to the location of the MINC2 file to
#' be read.
#' @return Returns a vector of mincSingleDim class
#' @seealso mincWriteVolume
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' testfile <- mincGetVolume("/tmp/rminctestdata/brain_cut_out.mnc")
#' }
#' @export
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
  attr(output, "sizes") <- sizes
  return(output)
}

#' @export
print.mincMultiDim <- function(x, ...) {
  cat("Multidimensional MINC volume\n")
  cat("Columns:      ", colnames(x), "\n")
  print(attr(x, "likeVolume"))
}

#' @export 
print.mincLogLikRatio <- function(x, ...) {
  cat("mincLogLikRatio output\n\n")
  mincLmerLists <- attr(x, "mincLmerLists")
  for (i in 1:length(mincLmerLists)) {
    cat("Model", i, ":\n")
    cat("  Formula:  ")
    print(mincLmerLists[[i]][[1]]$formula)
  }
  cat("\nMask used:", attr(x, "mask"), "\n")
  cat("Chi-squared Degrees of Freedom:", attr(x, "df"), "\n")
}

#' @export 
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
      

#' @export
print.vertexMultiDim <- function(x, ...) {
  cat("Multidimensional Vertstats file\n")
  cat("Columns:      ", colnames(x), "\n")
}
      
#' @export 
print.mincSingleDim <- function(x, ...) {
  cat("MINC volume\n")
  print(attr(x, "likeVolume"))
  print(attr(x, "filename"))
}

#' @export 
print.mincQvals <- function(x, ...) {
  print.mincMultiDim(x)
  cat("Degrees of Freedom:", paste(attr(x, "DF")), "\n")
  cat("FDR Thresholds:\n")
  print(attr(x, "thresholds"))
}


#' Volume Export
#' 
#' Writes a MINC volume to file
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
#' @param ... additional arguments to pass to methods
#' @details This function takes numeric data, usually the results computed
#' from one of the other mincFunctions, and writes it to file so that it
#' can be viewed or manipulated with the standard MINC tools
#' @return A list with the parameters of the minc volume written
#' @seealso mincWriteVolume,mincLm,mincFDR,mincMean,mincSd
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run a linear model relating the data in all voxels to Sex
#' vs <- mincLm(gf$jacobians_fixed_2 ~ Sex, gf)
#' # write the results to file
#' mincWriteVolume(vs, "Fstat.mnc", "F-statistic")
#' }
#' @export
mincWriteVolume <- function(buffer, ...) {
  UseMethod("mincWriteVolume")
}

#' @describeIn mincWriteVolume mincSingleDim
#' @export
mincWriteVolume.mincSingleDim <- function(buffer, output.filename, clobber = NULL, ...) {
  mincWriteVolume.default(buffer, output.filename, attr(buffer, "likeVolume"), clobber)
}

#' @describeIn mincWriteVolume mincMultiDim
#' @export
mincWriteVolume.mincMultiDim <- function(buffer, output.filename, column=1, 
																				like.filename = NULL, clobber = NULL, ...) {
  cat("Writing column", column, "to file", output.filename, "\n")
  if (is.null(like.filename)) {
    like.filename <- attr(buffer, "likeVolume")
  }
  if (is.na(file.info(as.character(like.filename))$size)) {
    stop(c("File ", like.filename, " cannot be found.\n"))
  }

  mincWriteVolume.default(buffer[,column], output.filename, like.filename, clobber)
}

#' @describeIn mincWriteVolume default
#' @export
mincWriteVolume.default <- function(buffer, output.filename, like.filename,
                                    clobber = NULL, ...) {
	
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

  if(length(which(is.nan(buffer))) != 0 || length(which(is.infinite(buffer))) != 0  || length(which(is.na(buffer))) != 0) {
	stop ("Cannot write volumes with inf,na or nans in them. Please remove the offending character")
   }


  output <- .C("write_minc2_volume",
               as.character(output.filename),
               as.character(like.filename),
               as.integer(start),
               as.integer(sizes),
               as.double(b.max),
               as.double(b.min),
               as.double(buffer), PACKAGE="RMINC")
  
  #Create starter history
  history <- 
    sprintf("%s>>> RMINC output volume, from a(n) %s buffer, like %s",
            Sys.time() %>% format("%a %b %d %H:%M:%S %Y"),
            class(buffer)[1],
            like.filename)
  
  #Initialize output volume's history to the starter history
  .Call("minc_overwrite_history", output.filename, history, nchar(history))
  
  return(invisible(NULL))
}
###########################################################################################


#' Dimension Sizes
#' 
#' Returns the dimension sizes of a MINC2 volume
#' 
#' 
#' @param filename The filename of the MINC2 volume whose dimension sizes are
#' to be returned.
#' @seealso minc.get.volume
#' @export
minc.dimensions.sizes <- function(filename) {
  sizes <- .C("get_volume_sizes",
              as.character(filename),
              sizes = integer(3), PACKAGE="RMINC")$sizes
  return(sizes)
}

#' Minc History
#' 
#' Retrieve or edit the history of a MINC Volume
#'
#' @param filename A path to a minc volume
#' @param new_history A new line to be added to the history
#' defaults to "[timestamp]>>> Written out by RMINC"
#' @return a character vector with one element per line of history
#' @name minc_history
NULL

#' describeIn minc_history retrieve
#' @export
minc.get.history <- 
  function(filename){
    history <- 
      .Call("get_minc_history",
            filename,
            20000) #Hardcode 20,000 character history limit - up for debate/fix
    
    unlist(strsplit(history, "\n"))
  }

#' describeIn minc_history append
#' @export
minc.append.history <-
  function(filename, 
           new_history = NULL){
    
    old_history <- minc.get.history(filename)
    
    if(is.null(new_history)){
      new_history <-
        sprintf("%s>>> Written out by RMINC",
                Sys.time() %>% format("%a %b %d %H:%M:%S %Y"))
    }
    
    new_history <- 
      paste0(c(old_history, new_history), collapse = "\n")
    
    .Call("minc_overwrite_history", filename, new_history, nchar(new_history))
    
    invisible(NULL)
  }


#' Get a hyperslab from a MINC2 file
#' 
#' Returns a 1D array by extracting a hyperslab of specified starts and counts
#' from a MINC2 volume.
#' 
#' This function allows for the extraction of an arbitrary contiguous chunk of
#' data from a MINC2 volume. The coordinates are voxel coordinates, given in
#' the volume dimension order.
#' 
#' @param filename The filename of the MINC2 volume from which to extract the
#' @param start A 3-dimensional array of voxel coordinate values to specify the
#' start of the hyperslab.
#' @param count A 3-dimensional array of voxel coordinate values to specify the
#' count of the hyperslab.
#' @return a numeric vector of size \code{prod(count)} containing the
#' hyperslab
#' @export
minc.get.hyperslab <- function(filename, start, count) {
#  total.size <- (count[1] - start[1]) * (count[2] - start[2]) *
#                (count[3] - start[3])
  total.size <- (count[1] * count[2] * count[3])
  buffer <- double(total.size)
  
  cat("Total size: "); cat(total.size); cat("\n")
  output <- .C("get_hyperslab",
               as.character(filename),
               as.integer(start),
               as.integer(count), hs=buffer, PACKAGE="RMINC")$hs
 
   return(buffer)
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

#' Voxel-wise ANOVA
#'  
#' Compute a sequential ANOVA at each voxel
#' @param formula The anova formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.
#' @param data The dataframe which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same length as the input files. ANOVA will only be computed
#' inside the mask.
#' @details This function computes a sequential ANOVA over a set of files.
#' @return Returns an array with the F-statistic for each model specified by formula with the following attributes: model – design matrix, filenames – 
#' 	minc file names input,dimensions,dimension names, stat-type: type of statistic used, df – degrees of freedom of each statistic. 
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run an ANOVA at each voxel
#' vs <- mincAnova(jacobians_fixed_2 ~ Sex, gf)
#' }
#' @export
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


#' @description Linear Model at Every Voxel
#' @name mincLm
#' @title Linear model at Every Voxel
#' @param formula The linear model formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.The RHS of the formula may contain one term with filenames. If so only the + operator may be used, and only two terms may appear on the RHS
#' @param data The data frame which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same length as the input files. The linear model will only be computed
#' inside the mask.
#' @param maskval the value in the mask used to select unmasked voxels, defaults to any positive intensity
#' from 1-99999999 internally expanded to .5 - 99999999.5. If a number is specified voxels with intensities 
#' within 0.5 of the chosen value are considered selected. 
#' @details This function computes a linear model at every voxel of a set of files. The function is a close cousin to lm, the key difference
#' being that the left-hand side of the formula specification takes a series of filenames for MINC files.
#' @return mincLm returns a mincMultiDim object which contains a series of columns corresponding to the terms in the linear model. The first
#' column is the F-statistic of the significance of the entire volume, the following columns contain the R-Squared term, the marginal t-statistics for each of the terms in the model along with their respective coefficients.
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # Compute a linear model at each voxel
#' vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
#' }
#' @export
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

  # Do Error Checking on mask. The mask may have a different dimension order than the input volumes
  # A call to minc.dimensions.sizes will read the dimension order according to the file order so 
  # it can be used to do the checking.
  if(!is.null(mask)) {
  	maskDim = minc.dimensions.sizes(mask)
  	volumeDim = minc.dimensions.sizes(parseLmOutput$data.matrix.left)
        # Case 1: Mask has one or more dimensions less than volume
	if(maskDim[1] != volumeDim[1] || maskDim[2] != volumeDim[2] || maskDim[3] != volumeDim[3]) {
		cat(paste("Mask Volume",as.character(maskDim[1]),as.character(maskDim[2]),as.character(maskDim[3]),
			  "Input Volume",as.character(volumeDim[1]),as.character(volumeDim[2]),as.character(volumeDim[3]),'\n'))
		stop("Mask Volume input volume dimension mismatch.") }
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
#' @export
pt2 <- function(q, df,log.p=FALSE) {
  2*pt(-abs(q), df, log.p=log.p)
}

#' Minc Masks
#' 
#' Either create a mask from a file name (by reading the volume)
#' or pass the vector along undisturbed.
#' 
#' @param mask Either the path to a mask file or a numeric vector representing the mask
#' @return a numeric mask vector
mincGetMask <- function(mask) {
  if (class(mask) == "character") {
    return(mincGetVolume(mask))
  }
  else {
    return(mask)
  }
}


#' Minc False Discovery Rates
#' 
#' Takes the output of a mincLm type run and computes the False Discovery Rate on the results.
#' @name mincFDR
#' @aliases vertexFDR anatFDR
#' @param buffer The results of a mincLm type run.
#' @param columns A vector of column names. By default the threshold will
#' be computed for all columns; with this argument the computation can
#' be limited to a subset.
#' @param mask Either a filename or a numeric vector representing a mask
#' only values inside the mask will be used to compute the
#' threshold.
#' @param df The degrees of freedom - normally this can be determined
#' from the input object.
#' @param statType This should be either a "t","F","u","chisq" or "tlmer" depending upon the
#' type of statistic being thresholded.
#' @param method The method used to compute the false discovery
#' rate. Options are "FDR" and "pFDR".
#' @param ... extra parameters to pass to methods
#' @details This function uses the \code{qvalue} package to compute the
#'  False Discovery Rate threshold for the results of a \link{mincLm}
#'  computation. The False Discovery Rate represents the percentage of
#'  results expected to be a false positive. Two implementations can be
#'  used as specified by the method argument. "FDR" uses the
#'  implementation in \code{p.adjust}, whereas "pFDR" is a version of the
#'  postivie False Discovery Rate as found in John Storey's \code{qvalue}
#'  package. The main interface functions are 
#'  \itemize{
#'  \item{mincFDR.mincMultiDim}{ The workhorse function, used to compute q-values
#'  and thresholds for sets of minc volumes}
#'  \item{mincFDR.logLikRatio}{ Similar to above, but calculates thresholds by parametric
#'  bootstrap when possible}
#'  \item{mincFDR.mincSingleDim}{ Used when \link{mincLm}-like results are written out and read back in
#'  either to the same or another R session. In this case it loses it's \code{minMultiDim} class
#'  and must be converted back}
#'  \item{vertexFDR}{ Used with results of a \link{vertexLm}-like command. Results are converted internally
#'  to resemble a mincMultiDim and processed as normal}
#'  } 
#' @return A object of type \code{mincQvals} with the same number of columns
#'  as the input (or the subset specified by the columns argument to
#'  mincFDR). Each column now contains the qvalues for each voxel. Areas
#'  outside the mask (if a mask was specified) will be represented by a
#'  value of 1. The result also has an attribute called "thresholds"
#'  which contains the 1, 5, 10, 15, and 20 percent false discovery rate
#'  thresholds.
#' @seealso mincWriteVolume,mincLm,mincWilcoxon or mincTtest 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
#' # compute the False Discovery Rate
#' qvals <- mincFDR(vs)
#' }
#' @export
mincFDR <- function(buffer, ...) {
  UseMethod("mincFDR")
}


#' @export
vertexFDR <- function(buffer, method="FDR", mask = mask) {
  mincFDR.mincMultiDim(buffer, columns=NULL, mask=NULL, df=NULL,
                       method=method)
}

#' @describeIn mincFDR mincSingleDim
#' @export
mincFDR.mincSingleDim <- function(buffer, df, mask=NULL, method="qvalue", ...) {
  if (is.null(df)) {
    stop("Error: need to specify the degrees of freedom")
  }
  if (length(df) == 1) {
    df <- c(1,df)
  }
  mincFDR.mincMultiDim(buffer, columns=1, mask=mask, df=df, method=method, ...)
}

#' mincFDRMask
#'
#' Returns either the specified mask, the mask associated with the buffer, 
#' or a vector of ones to be used as a mask.
#' 
#' @param mask a mask file or vector to be passed to mincGetMask
#' if left null, the buffer is checked for a mask attribute, if no mask
#' is found, a vector of ones is used
#' @param buffer a buffer describing a minc volume  
#' @return a numeric mask vector
#' @export
mincFDRMask <- function(mask = NULL, buffer) {
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
  return(mask)
}

#' @describeIn mincFDR mincLogLikRatio
#' @export
mincFDR.mincLogLikRatio <- function(buffer, mask=NULL, ...) {
  cat("Computing FDR for mincLogLikRatio\n")
  df <- attr(buffer, "df")
  mask <- mincFDRMask(mask, buffer)
  ncols <- ncol(buffer)
  
  # compute the thresholds at several sig levels
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  
  # check whether the parametric bootstrap was estimated on the model
  haveParametricBootstrap <- "parametricBootstrap" %in% names(attributes(buffer))
  if (haveParametricBootstrap) {
    # we'll keep the estimate, corrected q|pvals and their upper and lower 95% conf limits
    # note: same as parametricBootstrap, at the moment has a limit of a single column of
    # chi-squared values (i.e. the comparison of two input models)
    ncols <- 4
  }

  # compute p and q values
  # pvals and qvals size corresponds to voxels inside mask, whereas output
  # is the size of the buffer. (ncols is changed if parametricBootstrap present)
  pvals <- matrix(nrow=sum(mask>0.5), ncol=ncols)
  output <- matrix(1, nrow=nrow(buffer), ncol=ncols)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=ncols)
  qvals <- pvals
  

  for (i in 1:ncols) {
    # compute qvals through R's p.adjust function
    currentDF <- df[[1]]
    if (i == 1 | haveParametricBootstrap == F) {
      currentDF <- df[[i]]
      pvals[, i] <- pchisq(buffer[mask>0.5, i], currentDF, lower.tail=F)
    }
    if (haveParametricBootstrap & i==1) {
      # the linear model computed as part of the parametric bootstrap
      # note: as with the parametricBootstrap code, at the moment the assumption is
      # that there were only two models being compared, and thus a single column of
      # chisq values in the input 
      
      pvals[,2:4] <- predict(attr(buffer, "parametricBootstrapModel"),
                             newdata=data.frame(chisq=pvals[,i]),
                             interval="confidence")
    }
    qvals[, i] <- p.adjust(pvals[, i], "fdr")

    tfunc <- function(x) { qchisq(max(x), currentDF, lower.tail=FALSE) }
    thresholds[,i] <- mincFDRThresholdVector(pvals[,i], qvals[,i], tfunc, p.thresholds)
    output[mask>0.5, i] <- qvals[, i]
  }

  if (haveParametricBootstrap) {
    columnNames <- c("Chisq approx.", "Corrected", "Corrected lwr CI", "Corrected upr CI")
  }
  else {
    columnNames <- colnames(buffer)
  }
  
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columnNames
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columnNames
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()

  return(output)
}

#' a utility function to compute thresholds
#'
#' @param pvals a vector of pvalues
#' @param qvals a vector of corrected qvalues (such as returend by p.adjust)
#' @param thresholdFunc a function that returns the threshold given a vector of pvalues
#' @param p.thresholds the pvalues at which to compute the threshold
#'
#' The function should be the quantile function for the distribution being tested. For example,
#' for the chi squared distribution the function would be:
#' tfunc <- function(x) { qchisq(max(x), df[[i]], lower.tail=FALSE) }
mincFDRThresholdVector <- function(pvals, qvals, thresholdFunc=NULL,
                                   p.thresholds = c(0.01, 0.05, 0.10, 0.15, 0.20)) {
  
  thresholds <- vector("numeric", length=length(p.thresholds))
  for (j in 1:length(p.thresholds)) {
    # compute thresholds; to be honest, not quite sure what the NA checking is about
    subTholdPvalues <- pvals[qvals <= p.thresholds[j]]
    subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
    
    if ( length(subTholdPvaluesNumbers) >= 1 ) {
      thresholds[j] <- thresholdFunc(subTholdPvaluesNumbers)
        #qchisq(max(subTholdPvaluesNumbers), df[[i]], lower.tail=FALSE)
    }
    else { thresholds[j] <- NA }
  }
  return(thresholds)
}


#' @describeIn mincFDR mincLmer
#' @export
mincFDR.mincLmer <- function(buffer, mask=NULL, ...) {
  cat("In mincFDR.mincLmer\n")

  # if no DF set, exit with message
  df <- attr(buffer, "df")
  if (is.null(df)) {
    stop("No degrees of freedom for mincLmer object. Needs to be explicitly assigned with mincLmerEstimateDF (and read the documentation of that function to learn about the dragons that be living there!).")
  }
  else {
    warning("Here be dragons! Hypothesis testing with mixed effects models is challenging, since nobody quite knows how to correctly estimate denominator degrees of freedom.")
  }

  # get the mask
  mask <- mincFDRMask(mask, buffer)
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

#' @describeIn mincFDR mincMultiDim
#' @export
mincFDR.mincMultiDim <- function(buffer, columns=NULL, mask=NULL, df=NULL,
                                 method="FDR", statType=NULL) {
  
  if(method == "qvalue")  
    if(!requireNamespace("qvalue", quietly = TRUE))
      stop("The qvalue package must be installed for mincFDR to work")
    
  
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
      qobj <- qvalue::qvalue(pvals)
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


#' Minc Voxel Summary Functions
#' 
#' Compute the mean, standard deviation, sum, or variance at every voxel across a
#' a set of MINC volumes.
#' An optional grouping variable will split the computation by group
#' rather than performing it across all volumes as is the default.
#' @param filenames Filenames of the MINC volumes across which to create the
#' descriptive statistic.
#' @param grouping Optional grouping - contains same number of elements as
#' filenames; the results will then have the descriptive
#' statistic computed separately for each group.
#' @param mask A mask specifying which voxels are to be included in the
#' summary.
#' @param maskval the value in the mask used to select unmasked voxels, 
#' defaults to any positive intensity from 1-99999999 internally expanded to
#' .5 - 99999999.5. If a number is specified voxels with intensities 
#' within 0.5 of the chosen value are considered selected. 
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files. If a
#'          grouping factor was specified then the output will be a
#'          matrix consisiting of as many rows as there were voxels in
#'          the files, and as many columns as there were groups.
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mm <- mincMean(gf$jacobians_0.2) 
#' ms <- mincSd(gf$jacobians_0.2)
#' mv <- mincVar(gf$jacobians_0.2,gf$Strain) 
#' ms2 <- mincSum(gf$jacobians_0.2,gf$Strain)
#' }
#' @name mincSummaries 
NULL

#' @describeIn mincSummaries mean
#' @export
mincMean <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="mean", maskval=maskval)
  return(result)
}

#' @describeIn mincSummaries Variance
#' @export
mincVar <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  return(result)
}

#' @describeIn mincSummaries Sum
#' @export
mincSum <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="sum", maskval=maskval)
  return(result)
}

#' @describeIn mincSummaries Standard Deviation
#' @export
mincSd <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  result <- sqrt(result)
  return(result)
}

#' Minc T-test
#' 
#' Perform an unpaired,unequal variance t-test across a set of minc volumes
#' @param filenames Filenames of the MINC volumes across which to run the t-test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups with which to compare means
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files, with that voxel's t-statistic
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mtt <- mincTtest(gf$jacobians_0.2,gf$Strain)
#' } 
#' @export
mincTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  # the grouping for a t test should only contain 2 groups. Should 
  # also be converted to a factor if it's not.
  if( ! is.factor(grouping) ){
    grouping <- as.factor(grouping)
  }
  if(length(levels(grouping)) != 2 ){
    cat("\nThe groups (levels) in your data are:\n")
    cat(levels(grouping), "\n\n")
    stop("mincTtest can only be performed on data with 2 groups")
  }
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

#' @title Minc Paired T Test
#' @description Perform a paired t-test across a set of minc volumes
#' @param filenames Filenames of the MINC volumes across which to run the t-test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups with which to compare means. The two groups
#' must be the same length.
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files, with that voxel's t-statistic
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' gf = gf[1:20,]
#' mptt <- mincPairedTtest(gf$jacobians_0.2,gf$Strain)
#' }
#' @export
mincPairedTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  # here, similarly to the t-test, there should be 2 groups in the data. However
  # since this is a paired t-test (repeated measures, and it assumes element 1 from
  # group 1 is paired with element 1 from group 2 etc.), the groups need to have 
  # the same length.
  if( ! is.factor(grouping) ){
    grouping <- as.factor(grouping)
  }
  if(length(levels(grouping)) != 2 ){
    cat("\nThe groups (levels) in your data are:\n")
    cat(levels(grouping), "\n\n")
    stop("mincPairedTtest can only be performed on data with 2 groups")
  }
  lenght_group_1 <- sum(grouping == levels(grouping)[1])
  lenght_group_2 <- sum(grouping == levels(grouping)[2])
  if(lenght_group_1 != lenght_group_2) {
    cat("\nGroup 1 (",levels(grouping)[1], ") has ", lenght_group_1," elements.\n")
    cat("Group 2 (",levels(grouping)[2], ") has ", lenght_group_2," elements.\n\n")
    stop("The groups do not have the same length. mincPairedTtest expects 2 groups with equal length (paired data)")
  }
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


#' @title Minc Correlation
#' @description Perform a correlation between a set of minc volumes and a covariate.
#' @param filenames Filenames of the MINC volumes across which to run the correlation
#' @param grouping  Contains same number of elements as
#' filenames; contains values with which to correlate
#' @param mask A mask specifying which voxels are to be included in the correlation
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#' elements as there are voxels in the input files, with that voxel's correlation value (Pearson
#' correlation coefficient)
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mc <- mincCorrelation(gf$jacobians_0.2,gf$Weight)
#' }
#' @export 
mincCorrelation <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="correlation", maskval=maskval)
  return(result)
}

#' @title Minc Wilcoxon
#' @description Perform a Mann-Whitney U test between a set of minc volumes.
#' @param filenames Filenames of the MINC volumes across which to run the test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups.
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#' elements as there are voxels in the input files, with that voxel's U value (lower one).
#' The number of observations in each sample is also saved as an attribute(m and n) so the result
#' can be passed into mincFDR.
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mw <- mincWilcoxon(gf$jacobians_0.2,gf $Strain)
#' }
#' @export
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

#' @describeIn mincApply parallel
#' @export
pMincApply <- function(filenames, function.string,
                       mask=NULL, workers=4, tinyMask=FALSE, 
                       method="snowfall",global="",packages="", 
                       modules="",vmem="8",walltime="01:00:00") {
  
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
  
  # Saving to /tmp does not always work...
  maskFilename <- paste("pmincApplyTmpMask-", Sys.getpid(), ".mnc", sep="")
  
  #If the current working directory isn't writeable, 
  #write to a tempdir instead
  if(file.access(getwd(), 2) != 0) maskFilename <- file.path(tempdir(), maskFilename)
  
  mincWriteVolume(maskV, 
                  maskFilename, 
                  clobber=TRUE) 
  
  
  # create the packageList that will be used for the snowfall and sge options
  # if packages contains multiple libraries, the test (packages == "") 
  # will return as many TRUE/FALSE as the length of the vector. So to test
  # for "", first test that the length of the packages vector is 1
  if(length(packages) < 2) {
    if(packages == "") {
      packageList = c("RMINC")
    }
    else {
      packageList = c(packages,"RMINC")
    }
  }
  else {
    packageList = c(packages,"RMINC")
  }
  
  pout <- list()
  
  if (method == "local") {
    stop("Lovely code ... that generates inconsistent results because something somewhere is not thread safe ...")
    
    if(!(requireNamespace("doMC", quietly = TRUE) & requireNamespace("foreach", quietly = TRUE))) 
      stop("One or both of doMC and foreach is missing, please install these packages")
    
    registerDoMC(workers)
    
    # run the job spread across each core
    pout <- foreach(i=1:workers) %dopar% { mincApply(filenames, function.string,
                                                   mask=maskFilename, maskval=i) }
    #cat("length: ", length(pout), "\n")
  }

  # The pbs options use mpirun and snow to parallelize mincApply over multiple cores.
  # It is currently configured to send all the jobs to one node, and parallelize over the
  # cores at that node. Also note the amount of virtual memory requested is set at 8g.
 
  # It operates as follows:
  # 1) Save global variables to disk
  # 2) Write out a .R file that will execute mpi operations once submitted to the cluster
  # 3) Write out a .sh file that will be submitted to the cluster
  # 4) Submit the .sh file to the cluster
  # 5) Wait for the jobs to finish
  # 6) Read the output from disk

  else if (method == "pbs") {
 
   if(is.null(getOption("MAX_NODES"))) {
		 ppn = 8
                 #maximize node usage
	 	 nodes = ceiling(workers/ppn)
		 workers = nodes*ppn-1
   }
   else {
		 ppn = getOption("MAX_NODES")
		 nodes = ceiling(workers/ppn)
    }

    # 1) Save variables which will be referenced on the cluster to disk (including user specified global variables)
    rCommand = 'save(\'maskV\',\'filenames\',\'workers\',\'REDUCE\',\'function.string\',\'maskFilename\','
    for (nVar in 1:length(global)) {
	if(global[nVar] != "") {
		rCommand = paste(rCommand,'\'',global[nVar],'\',',sep="")
        }
    }
    rCommand = paste(rCommand,'file=\'mpi-rminc-var\')',sep="")
    
    eval(parse(text=rCommand))

    # 2) Write out an R file to disk. This file will be executed via mpirun on the cluster
    fileConn <- file("mpi-rminc.R",open='w')

    # Load R Packages
    

    # snow is used to coordinate operations, but Rmpi could be used for greater control
    packageList = c(packageList,"snow")

    for (nPackage in 1:length(packageList)) {     
    	writeLines(paste("library( ",packageList[nPackage],")",sep = ""),fileConn)
    }

    # Load the variables we saved in step 1
    writeLines("load(\"mpi-rminc-var\")",fileConn)

    # Create the snow cluster
    writeLines("cl <- makeCluster(workers, type = \"MPI\")",fileConn) 

    # Create a wrapper function that we can easily run with clusterApply
    writeLines("wrapper <- function(i) { return(mincApply(filenames, function.string, mask = maskFilename, maskval = i, reduce = REDUCE))}",fileConn)

    # Export all neccessary variables to each slave
    writeLines("clusterExport(cl,c('filenames','REDUCE','function.string','maskFilename','maskV'))",fileConn)

    # Main mpi exection
    writeLines("clusterOut <- clusterApply(cl,1:workers,wrapper)",fileConn)

    # At this point we are done.
    writeLines("stopCluster(cl)",fileConn)

    # Test data at one voxel to determine how many ouytputs
    writeLines(" x <- mincGetVoxel(filenames, 0,0,0)",fileConn)
    writeLines("test <- eval(function.string)",fileConn)
    writeLines("if (length(test) > 1) {",fileConn)
    writeLines("output <- matrix(0, nrow=length(maskV), ncol=length(test))",fileConn)
    writeLines("class(output) <- class(clusterOut[[1]])",fileConn)
    writeLines("attr(output, \"likeVolume\") <- attr(clusterOut[[1]], \"likeVolume\")",fileConn)
    writeLines("} else {",fileConn)
    writeLines("output <- maskV",fileConn)
    writeLines("}",fileConn)
    writeLines("for(i in 1:workers) {",fileConn)
    writeLines("if (length(test)>1) {",fileConn)
    writeLines("if(REDUCE == TRUE)",fileConn)	
    writeLines("output[maskV == i,] <- clusterOut[[i]]",fileConn)
    writeLines("else",fileConn)
    writeLines("output[maskV == i,] <- clusterOut[[i]][maskV == i, ]",fileConn)
    writeLines("}",fileConn)
    writeLines("else {",fileConn)
    writeLines("output[maskV==i] <- clusterOut[[i]]",fileConn)
    writeLines("}",fileConn)
    writeLines("}",fileConn)

    # Write the output (R data) to disk, so the initiating R session can access the data
    writeLines("save('output',file='mpi-rminc-out')",fileConn)

    close(fileConn)

    # 3) Write out a .sh file which will be what is submitted to the cluster
    qfileConn <- file("q-mpi-rminc.sh",open='w')
    writeLines("#!/bin/bash -x",qfileConn)
   
    # Errors and Output are written to the current directory, but this could be set to a user defined spot.
    writeLines("#PBS -e ./",qfileConn)
    writeLines("#PBS -o ./",qfileConn)
    writeLines("#PBS -N pMincApply",qfileConn) 

    # Allocate nodes and ppn
    writeLines(paste("#PBS -l nodes=",as.character(nodes),":ppn=",as.character(ppn),sep=""),qfileConn)

    # Allocate walltime and vmem
    writeLines(paste("#PBS -l vmem=",vmem,"g,walltime=",walltime,sep=""),qfileConn)

    # Load modules
    for (nModule in 1:length(modules)) {     
        if(modules[nModule] != "") {
    		writeLines(paste("module load ",modules[nModule],sep = ""),qfileConn)
       }
    }

    # Define temp directory
    writeLines(paste("export TMPDIR=",getOption("TMPDIR"),sep=""),qfileConn)

   # For Testing
    writeLines(paste("cp -R ~/Software/RMINC/rminctestdata /tmp/",sep=""),qfileConn)

    # Move to working directory
    writeLines(paste("cd ",getOption("WORKDIR"),sep=""),qfileConn)

    # Neccessary to initiate mpi operations.
    writeLines("mpirun -np 1 R CMD BATCH mpi-rminc.R",qfileConn)

    close(qfileConn)

    # 4) Submit the job
    result <- system("qsub q-mpi-rminc.sh",intern=TRUE)
    ptm = proc.time()
    
    # 5) Wait for Job to finish
    status <- system(paste("qstat ",result," | grep C"),intern=TRUE) 
    # 1 Indicates not completed
    while(length(status) == 0) {
       flush.console()
       status <- system(paste("qstat ",result," | grep R"),intern=TRUE)
       runTime = proc.time()-ptm
       if(length(status) == 0) {
	  cat(paste("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bQueued:  ",sprintf("%3.3f seconds",(runTime[3]))))
       } else {
	  cat(paste("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRunning: ",sprintf("%3.3f seconds",(runTime[3]))))
       }
       Sys.sleep(2)
       status <- system(paste("qstat ",result," | grep C"),intern=TRUE)
      
    }
    print(paste("Completed: ",as.character(runTime[3])," Seconds"))
    # 6) Read output from disk
    load("mpi-rminc-out")
    
    # Clean up intermediate files
    system(paste("rm",maskFilename))
    system("rm q-mpi-rminc.sh")
    system("rm mpi-rminc.R")
    system("rm mpi-rminc-var")
    system("rm mpi-rminc-out")

    return(output)

}

  else if (method == "sge") {
    
    if(!requireNamespace("Rsge", quietly = TRUE)) 
      stop("The Rsge package is require to run code with SGE parallelism")

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

    # Need to use double quotes, because both sge.submit and mincApply try to evalute the function
    function.string = enquote(function.string)
    
    l1 <- list(length=workers)

    if(global == "")  {
      globallist = c(sub("\\(([A-Z]|[a-z])\\)","",function.string))
    }
    else {
      globallist = c(global,sub("\\(([A-Z]|[a-z])\\)","",function.string))
    }

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
    
    if(!requireNamespace("snowfall", quietly = TRUE))
      stop("The snowfall package is required to run code with snowfall parallelism")
    
    sfInit(parallel=TRUE, cpus=workers)

    for (nPackage in 1:length(packageList)) {
      sfLibrary(packageList[nPackage],character.only=TRUE) 
    }
		
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


#' @description Can execute any R function at every voxel for a set of MINC files
#' @name mincApply
#' @title Apply arbitrary R function at every voxel
#' @param filenames The MINC files over which to apply the function. Have to be
#' the same sampling.
#' @param function.string The function which to apply. Can only take a single
#' argument, which has to be 'x'. See details and examples.
#' @param mask The filename of a mask - function will only be evaluated
#' inside the mask.
#' @param maskval An integer specifying the value inside the mask where to
#' apply the function. If left blank (the default) then anything
#' above 0.5 will be considered inside the mask. This argument
#' only works for mincApply, not pMincApply.
#' @param workers The number of processes to split the function into. Only
#' works for pMincApply, not for mincApply.
#' @param method The method to be used for parallization ("snowfall" ,"sge","scinet","hpf")
#' "hpf" : The SickKids Cluster --> Currently 8g of virtual memory is requested. Additionally
#' only one node (32 workers) can be requested
#' "scinet" : scinet cluster --> scinet allows mutliple nodes to be requested but currently the limit is only one node
#' (max 8 workers). The walltime is also set to 1:00:00 and virtual memory is also 8g
#' @param tinyMask Only for pMincApply; if set to some numeric value it computes
#' the function over that fraction of the mask. Useful for
#' debugging purposes only (i.e. if you want to test out whether
#' a new function works across the cluster.)
#' @param modules modules to be loaded on the cluster nodes
#' @param vmem Amount of virtual memory to request
#' @param walltime Amount of walltime to request
#' @param reduce Whether or not to collapse the results into a \code{mincMultiDim} object
#' @param global Global variables to be exported
#' @param packages R Packages to be exported
#' @details mincApply allows one to execute any R function at every voxel of a
#' set of files. There are two variants: mincApply, which works
#' inside the current R session, and pMincApply, which uses the
#' snowfall and Rsge packages to split the execution of
#' the function across multiple cores/processors on the same machine
#' or across a cluster.
#' Unless the function to be applied takes a single argument a
#' wrapper function has to be written. This is illustrated in the
#' examples; briefly, the wrapper function takes a single argument,
#' called 'x', which will take on the voxel values at every voxel.
#' The function has to be turned into a string; the quote function
#' can be very handy. The output of this function also has to be a
#' simple vector or scalar.
#' Note that interpreted R can be very slow. Mindnumbingly slow. It
#' therefore pays to write optimal functions or, whenever available,
#' use the optimized equivalents. In other words and to give one
#' example, use mincLm rather than applying lm, and if lm has to
#' really be applied, try to use lm.fit rather than plain lm.
#' When using the pbs method, one can also set the options --> TMPDIR,MAX_NODES and WORKDIR
#' @return out The output is a matrix with the same number of rows a the
#' file sizes and the same number of columns as output by the
#' function that was applied. Cast into one of the MINC classes
#' to make writing it out with mincWriteVolume easier.

#' @seealso mincWriteVolume,mincMean,mincSd,mincLm,mincAnova
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' ma <- mincApply(gf$jacobians_fixed_2,quote(mean(x))); 
#' mincWriteVolume(ma, "means.mnc")
#' ### run the whole thing in parallel on the local machine
#' testFunc <- function (x) { return(c(1,2))}
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    workers = 4,
#'                    global = c('gf','testFunc'))
#' mincWriteVolume(pout, "pmincOut.mnc")
#' 
#' ### pbs example 1 (hpf)
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' testFunc <- function (x) { return(c(1,2))}
#' options(TMPDIR="/localhd/$PBS_ID")
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    workers = 4,
#'                    method="pbs",
#'                    global = c('gf','testFunc'))
#'
#' ### pbs example 2 (scinet)
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' testFunc <- function (x) { return(c(1,2))}
#' options(WORKDIR="$SCRATCH")
#' options(MAX_NODES=8)
#' options(TMP_DIR="/tmp")
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    modules=c("intel","openmpi","R/3.1.1"),
#'                    workers = 4,
#'                    method="pbs",
#'                    global = c('gf','testFunc'))
#'                    
#' #NOTE 1: On SCINET jobs are limited to 32*48 hours of cpu time
#' #NOTE 2: On SCINET the Rmpi library must be compiled for use with R 3.1.1
#' }
#' @export
mincApply <- 
  function(filenames, function.string, mask=NULL, maskval=NULL, reduce=FALSE) {
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

#' @export
vertexTable <- function(filenames) {
  nSubjects <- length(filenames)
  nvertices <- nrow(read.table(filenames[1]))
  output <- matrix(nrow=nvertices, ncol=nSubjects)
  for (i in 1:nSubjects) {
    output[,i] <- as.matrix(read.table(filenames[i]))
  }
  return(output)
}


#' Performs ANOVA on each vertex point specified 
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param subset rows to be used, by default all are used
#' @return Returns an array with the F-statistic for each model specified by formula with the following attributes: model – design matrix, filenames – 
#' 	vertex file names input, stat-type: type of statistic used,dimensions,dimension names, and df – degrees of freedom of each statistic. 
#' @seealso mincAnova,anatAnova 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' result = vertexAnova(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf)
#' }
#' @export 
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


#' Calculates statistics and coefficients for linear model of specified vertex files
#' @param formula a model formula. The RHS of the formula may contain one term with filenames. If
#' so only the + operator may be used, and only two terms may appear on the RHS
#' @param data a data.frame containing variables in formula 
#' @param subset rows to be used, by default all are used
#' @return Returns an object containing the R-Squared value,beta coefficients, F 
#' and t statistcs that can be passed directly into vertexFDR.
#' @seealso mincLm,anatLm,vertexFDR 
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' result = vertexLm(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf) 
#' }
#' @export
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

###########################################################################################
#' @description This function is used to compute an arbitrary function of every region in a set of vertex files.
#' @name vertexApply
#' @title Apply function over vertex Files
#' @param filenames vertex file names
#' @param function.string The function which to apply. Can only take a single
#' argument, which has to be 'x'.
#' @return  out: The output will be a single vector containing as many
#'          elements as there are vertices in the input variable. 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- vertexApply(gf$CIVETFILES$nativeRMStlink20mmleft,quote(mean(x)))
#' }
#' @export
vertexApply <- function(filenames,function.string) 
{
  # Load the data
  vertexData = vertexTable(filenames)

  # In order to maintain the same interface as mincApply, the (x) part needs to be stripped
  function.string = gsub('(x)','', function.string, fixed = TRUE)

  # The apply part (transpose to match output of mincApply)
  results <- t(apply(vertexData,1,function.string[1]))

  attr(results, "likeFile") <- filenames[1]
  
  # run the garbage collector...
  gcout <- gc()
  
  return(results)
}


#' Create descriptive statistics across a series of vertex files
#' 
#' This function is used to compute the mean, standard deviation, sum, or
#' variance of every vertex in a set of vertex files.
#' 
#' 
#' @param filenames Filenames of the vertex volumes across which to create the
#' descriptive statistic.
#' @return \item{out}{The output will be a single vector containing as many
#' elements as there are vertices in the input files.}
#' @seealso vertexLm
#' @examples
#' 
#' \dontrun{
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # compute the mean at every voxel of all files.
#' means <- vertexMean(gf$filenames)
#' }
#' @name vertexSummaries
NULL

#' @describeIn vertexSummaries mean
#' @export
vertexMean <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(rowMeans(vertexData))

} 

#' @describeIn vertexSummaries sum
#' @export
vertexSum <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(rowSums(vertexData))

} 

#' @describeIn vertexSummaries var
#' @export
vertexVar <- function(filenames) 
{
	vertexData = vertexTable(filenames)
	return(apply(vertexData,1,var))

} 

#' @describeIn vertexSummaries standard deviation
#' @export
vertexSd<- function(filenames) 
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
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' writeVertex(gf$nativeRMStlink20mm,"~/RMStlink20mm.txt",FALSE,NULL,NULL)
#' }
#' @export
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
        if (is.data.frame(gf)) {
            write("<matrix>", file = filename, append = TRUE)
            write.table(gf, file = filename, append = TRUE, 
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
#' @seealso minc.apply, minc.write.volume, minc.model
#' @examples
#' 
#' \dontrun{
#' # get a file that could be used by glim_image
#' gf <- as.data.frame(read.table("filename.glim"))
#' 
#' # get a t-test at every voxel
#' t.stats <- minc.model(gf$V1, gf$V2, "t-test")
#' 
#' # and write to file
#' minc.write.volume("t-test.mnc", gf$V1[1], t.stats)
#' 
#' # create a pretty picture to include in the next Nature or Science
#' # article
#' minc.ray.trace("t-test.mnc", output="pretty.rgb", threshold=c(2.5,6),
#'                background="anatomy_image.mnc")
#' }
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
#' @details mincLmer provides an interface to running linear mixed effects models at every
#' voxel. Unlike standard mincLm, however, testing hypotheses in linear mixed effects models
#' is more difficult, since the denominator degrees of freedom are more difficult to
#' determine. RMINC provides two alternatives: (1) estimating degrees of freedom using the
#' \code{\link{mincLmerEstimateDF}} function, and (2) comparing two separate models using
#' \code{\link{mincLogLikRatio}} (which in turn can be corrected using
#' \code{\link{mincLogLikRatioParametricBootstrap}}). For the most likely models - longitudinal
#' models with a separate intercept or separate intercept and slope per subject - both of these
#' approximations are likely correct. Be careful in using these approximations if
#' using more complicated random effects structures.
#'
#' @seealso \code{\link{lmer}} for description of lmer and lmer formulas; \code{\link{mincLm}}
#' @examples
#' \dontrun{
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc")
#' mincWriteVolume(vs, "age-term.mnc", "tvalue-age")
#' # run in parallel with multiple processors on the local machine
#' vs <- mincLmer(filenames ~ age + sex + (age|id), 
#'                data=gf, 
#'                mask="mask.mnc", 
#'                parallel=c("snowfall", 4))
#' # run in parallel with multiple processors over the sge batch queueing system
#' vs <- mincLmer(filenames ~ age + sex + (age|id), 
#'                data=gf, 
#'                mask="mask.mnc", 
#'                parallel=c("sge", 4))
#' # estimate degrees of freedom
#' vs <- mincLmerEstimateDF(vs)
#' # correct for multiple comparisons using the False Discovery Rate
#' (qvals <- mincFDR(vs))
#' # generate another model with a more complex curve for the age term
#' library(splines)
#' vs2 <- mincLmer(filenames ~ ns(age,2) + sex + (age|id), data=gf, mask="mask.mnc")
#' # see if that more complex age term was worth it
#' modelCompare <- mincLogLikRatio(vs, vs2)
#' mincFDR(modelCompare)
#' # see if there was any bias in those p-value estimates (takes a few minutes)
#' modelCompare <- mincLogLikRatioParametricBootstrap(modelCompare)
#' mincFDR(modelCompare)
#' }
#' @export
mincLmer <- function(formula, data, mask=NULL, parallel=NULL,
                     REML=TRUE, control=lmerControl(), start=NULL, verbose=0L) {

  #Try and load lme4
#   result = tryCatch({
# 	library(lme4)
#   }, error = function(e) {
# 	stop("Could not find lme4. Please install this package.")
#   })

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
#' of the model (see here \url{http://glmm.wikidot.com/faq#df}). mincLmer by
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
#' @export
mincLmerEstimateDF <- function(model) {
  # set the DF based on the Satterthwaite approximation
  # load lmerTest library if not loaded; lmerTest takes over some lmer functions, so unload if
  # it wasn't loaded in the first place
  # lmerTestLoaded <- "package:lmerTest" %in% search()
  # if (lmerTestLoaded == FALSE) {
  #   library(lmerTest)
  # }

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
  # if (lmerTestLoaded == FALSE) {
  #   detach("package:lmerTest")
  # }
  return(model)
}
# the actual optimization of the mixed effects models; everything that has to be recomputed
# for every voxel goes here. Works on x (each voxel is assigned x during the loop), and
# assumes that all the other info is in a variable called mincLmerList in the global
# environment. This last part is a hack to get around the lack of multiple function arguments
# for mincApply and friends.
#' @export
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
#' @export
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
#' @export
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

#' @export 
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
#' @param ... Two or more mincLmer objects
#' @return the voxel wise log likelihood test. Will have a number of columns corresponding to the
#' number of inputs -1. Note that it resorts the inputs from lowest to highest degrees of freedom
#'
#' @seealso \code{\link{lmer}} and \code{\link{mincLmer}} for description of lmer and mincLmer.
#' \code{\link{mincFDR}} for using the False Discovery Rate to correct for multiple comparisons,
#' and \code{\link{mincWriteVolume}} for outputting the values to MINC files.
#' @examples
#' \dontrun{
#' m1 <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc", REML=F)
#' m2 <- mincLmer(filenames ~ age + I(age^2) + sex + (age|id), 
#'                data=gf, mask="mask.mnc", REML=F)
#' m3 <- mincLmer(filenames ~ age + I(age^2) + I(age^3) + sex + (age|id), 
#'                data=gf, mask="mask.mnc", REML=F)
#' llr <- mincLogLikRatio(m1, m2, m3)
#' mincFDR(llr)
#' mincWriteVolume(llr, "m2vsm3.mnc", "m3")
#' }
#' @export
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
#' @return a matrix containing the chi-square p-values and the bootstrapped p-values
#' @export
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
  # create a linear model of the relation between the chi square assumption and
  # the parametric bootstrap (note the lack of intercept - when p is exactly 0
  # it should be 0 in both cases
  lmodel <- lm(parametricBootstrap ~ chisq -1, data=data.frame(out))
  # make the parametric bootstrap estimates an attribute of the logLik output
  attr(logLikOutput, "parametricBootstrap") <- out
  # make the model estimate an attribute as well
  attr(logLikOutput, "parametricBootstrapModel") <- lmodel
  return(logLikOutput)
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
#' @export
#' @examples
#' \dontrun{
#' index <- mincVectorToVoxelCoordinates("filename.mnc", 345322)
#' voxel <- mincGetVoxel(gf$filenames, index)
#' }
mincVectorToVoxelCoordinates <- function(volumeFileName, vectorCoord) {
  sizes <- minc.dimensions.sizes(volumeFileName)
  # the fun off by one bit
  index <- vectorCoord-1
  i3 <- index %% sizes[3]
  i2 <- (index / sizes[3]) %% sizes[2]
  i1 <- ((index / sizes[3]) / sizes[2]) %% sizes[1]
  return(floor(c(i1, i2, i3)))
}

#' selects a few random indices from a volume
#'
#' Given a filename, select a few random indices using the uniform distribution
#' from voxels that have a value of 1 (i.e. from a mask volume)
#'
#' @param volumeFileName the filename for a MINC volume
#' @param nvoxels the number of voxels to select
#' @param convert whether to convert to MINC voxel space (default) or keep in index space
#' @param ... additional arguments  
#' @return A vector of length \code{nvoxels} containing selected voxel indices or if convert is true
#' a matrix containing the x-y-z coordinates of the selected voxels. 
#' @export
mincSelectRandomVoxels <- function(volumeFileName, nvoxels=50, convert=TRUE, ...) {
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

#' @title Run Testbed
#' @description Run the test bed to ensure all RMINC functions
#' work on your system
#' @param verboseTest
#' Whether or not to verbosely print test output, default is
#' to print simplified results
#' @param purgeData whether to remove downloaded test files
#' in /tmp/rminctestdata
#' @param ... additional parameter for \link[testthat]{test_dir}
#' @return invisibly return the test results
#' @export
runRMINCTestbed <- function(..., verboseTest = FALSE, purgeData = TRUE) {
  
	# if(!require(testthat)){
	#   stop("Sorry, you need to install testthat to run the testbed")
	# }

  options(verbose = verboseTest)
  # Make sure environment is clear
  #rm(list=ls())
  
  if(!file.exists("/tmp/rminctestdata/")){
    system('mkdir /tmp/rminctestdata')
  }
  # Download Tarball from Wiki
  if(!file.exists("/tmp/rminctestdata/rminctestdata.tar.gz")){
    system("wget -O /tmp/rminctestdata/rminctestdata.tar.gz --no-check-certificate https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz")
  }
  # Untar
  system('tar -xf /tmp/rminctestdata/rminctestdata.tar.gz -C /tmp/')

  # Run Tests
  rmincPath = find.package("RMINC")
  cat("\n\nRunning tests in: ", paste(rmincPath,"/","tests/",sep=""), "\n\n\n")
  testReport <- testthat::test_dir(paste(rmincPath,"/","tests/",sep=""), ...)
  
  cat("\n*********************************************\n")
  cat("The RMINC test bed finished running all tests\n")
  cat("*********************************************\n\n\n")
  
  if(purgeData){
    cat("Removing temporary directory /tmp/rminctestdata\n")
    system('rm -fr /tmp/rminctestdata')
  }

  return(invisible(testReport))
}

#' Perform Arbitrary calculations on a collection of mincVolumes
#' 
#' An RCPP variant of \link{mincApply}, the primary advantage being
#' that functions of an arbitrary number of arguments can be passed
#' to mincApplyRCPP.
#' @param filenames The name of the files to apply over
#' @param fun the function to apply
#' @param ... additional parameters to fun
#' @param mask a numeric mask vector
#' @param maskval An integer specifying the value inside the mask where to
#' apply the function. If left blank (the default) then anything
#' above 0.5 will be considered inside the mask. This argument
#' only works for mincApply, not pMincApply.
#' @param filter_masked Whether or not to remove the masked values
#' from the resultant object
#' @param collate A function to (potentially) collapse the result list
#' examples include link{unlist} and \link{simplify2array}, defaulting
#' to \link{identify} which returns the unaltered list.
#' @return a list of results subject the the collate function
#' @export
mincApplyRCPP <- 
  function(filenames, 
           fun, 
           ..., 
           mask = NULL, 
           maskval = NULL,
           filter_masked = FALSE,
           collate = identity){
  
  apply_fun <- 
    function(x, extra_arguments)
      do.call(match.fun(fun), c(list(x), extra_arguments))
  
  args <- list(...)
  filenames <- as.character(filenames)
  
  if(is.null(maskval)){
    minmask <- 1
    maxmask <- 99999999
  } else {
    minmask <- maxmask <- maskval
  }
  
  if(is.null(mask)){
    use_mask <- FALSE
    mask = ""
  } else {
    use_mask <- TRUE
  }
  
  masked_value <- options()$RMINC_MASKED_VALUE
  
  results <-
    .Call("RMINC_rcpp_minc_apply",
        filenames,
        use_mask = use_mask,
        mask = mask,
        mask_lower_val = minmask,
        mask_upper_val = maxmask,
        value_for_mask = masked_value,
        filter_masked = filter_masked,
        fun = apply_fun,
        args = args,
        PACKAGE = "RMINC")
  
  # run the garbage collector...
  gcout <- gc()
  
  collation_function <- match.fun(collate)
  results <- collation_function(results)
  attr(results, "likeVolume") <- filenames[1]
  attr(results, "filenames") <- filenames
  
  return(results)
}

#' Download Example Data
#' 
#' Download the example data needed to run our examples in your /tmp directory
#' The data can be downloaded manually from 
#' \url{https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz}
#' @export
getRMINCTestData <- function() {

  system('mkdir /tmp/rminctestdata')

  # Download Tarball from Wiki
  system("wget -O /tmp/rminctestdata/rminctestdata.tar.gz --no-check-certificate https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz")

  # Untar
  system('tar -xf /tmp/rminctestdata/rminctestdata.tar.gz -C /tmp/')
  
}

# Run function with/without output silenced; used in test bed
#' @export
verboseRun <- function(expr,verbose,env = parent.frame()) {
	
	env$expr <- expr
	
	if(!verbose) {
	  sink("/dev/null")
	  on.exit(sink())
	}
	
	output = with(env,eval(parse(text=expr)))
	return(invisible(output))
}
