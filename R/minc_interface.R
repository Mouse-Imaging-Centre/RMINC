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
#' @export
print.RMINC_MASKED_VALUE <-
  function(x, ...)
    print(as.symbol("masked"))

# Attributes for minc result objects
known_minc_attributes <-
  c("filename", "filenames", "likeVolume", "mask", "mincLmerList", "mincLmerLists", "sizes",
    "model", "stat-type", "df", "dimnames")

#' Get and Set Minc Specific Attributes
#' 
#' Manage auxillary information contained in RMINC objects. Currently recognizes:
#' \itemize{
#' \item{filenames: filenames of minc files used to create an object}
#' \item{likeVolume: a minc file to use a structural template when writing out
#' data}
#' \item{mask: a mask file associated with the object}
#' \item{mincLmerList: extra information used in fitting \link{mincLmer}s}
#' \item{mincLmerLists: collections of mincLmerLists used to compare \link{mincLmer} fits}
#' }
#' @param minc_object the RMINC object of interest typically a \code{mincMultiDim} or
#' a related class
#' @param updated_attrs A named list containing the new attributes to be replaced
#' @return \code{mincAttributes} returns an attribute list, \code{setMincAttributes} returns 
#' the updated object
#' @export
mincAttributes <-
  function(minc_object){
    minc_attributes <- attributes(minc_object)
    minc_attributes <- 
      minc_attributes[names(minc_attributes) %in% known_minc_attributes]
    
    return(minc_attributes)
  }

#' @describeIn mincAttributes setter
#' @export
setMincAttributes <-
  function(minc_object, updated_attrs){
    
    if(is.null(updated_attrs) || length(updated_attrs) == 0)
      return(minc_object)
    
    if(any(! names(updated_attrs) %in% known_minc_attributes)) 
      stop(sprintf("New attributes must be a known minc object attribute (%s)",
                   paste0(known_minc_attributes, collapse = ", ")))
    
    all_attrs <- attributes(minc_object)
    
    if(is.null(all_attrs)){
      all_attrs <- environment()
    } else  {
      all_attrs <- as.environment(all_attrs)
    }
    
    list2env(updated_attrs, all_attrs)
    attributes(minc_object) <- as.list(all_attrs)
    
    return(minc_object)
  }

#' Collate Minc
#' 
#' Helper function to collate the results of a \link{mincApplyRCPP} family 
#' (\link{pMincApply}, \link{mcMincApply}, and \link{qMincApply}) function
#' @param result_list The mincApply results to collate.
#' @return a matrix like object of class \code{mincSingleDim}, code{mincMultiDim},
#' or code{mincList} depending on the dimensions of the input object
#'@export
simplify2minc <- function(result_list){
  
  ## Deal with masking, find the first non-masked value, NA it out
  ## insert it back in to standardize element length as much as possible
  lgl_missing <- vapply(result_list, 
                        function(res){ inherits(res, "RMINC_MASKED_VALUE") | all(is.na(res)) }, 
                        logical(1))
  first_element <- result_list[[min(which(!lgl_missing))]]      
  na_value <- first_element                                
  na_value[] <- getOption("RMINC_MASKED_VALUE")            #set all its elements to masked
  result_list[lgl_missing] <- list(na_value)                #replace in result list
  
  ## Determine the correct reduction technique and apply it
  if(first_element %>% is.atomic && length(first_element) == 1){
    simplified_results <- unlist(result_list)
  } else if(first_element %>% is.atomic && length(first_element) > 1){
    simplified_results <- t(simplify2array(result_list))
  } else if(first_element %>% is.data.frame) {
    simplified_results <- bind_rows(result_list)
  } else {
    simplified_results <- result_list
  }
  
  ## If any known_minc_attributes were already defined for the object
  ## ensure they carry over to the simplified object
  result_attributes <- mincAttributes(result_list)
  if(!is.null(result_attributes))
    simplified_results <- setMincAttributes(simplified_results, result_list)
  
  ## Reclass the object to the appropriate RMINC class
  if(is.null(ncol(simplified_results)) && is.list(simplified_results)){
    class(simplified_results) <- "mincList"
  } else if(is.null(ncol(simplified_results))){
    class(simplified_results) <- c("mincSingleDim", class(simplified_results))
  } else {
    class(simplified_results) <- c("mincMultiDim", class(simplified_results))
  }

  return(simplified_results)
}

#' Check for RMINC objects
#' 
#' Checks if the current object is one produced by RMINC
#' (currently \code{mincList}, \code{mincSingleDim}, \code{mincMultiDim})
#' @param x the object of interest
#' @return logical whether or not the object is an RMINC object
#' @export
is.minc <- function(x){
  inherits(x, c("mincList", "mincSingleDim", "mincMultiDim"))
}

#' Coerce to RMINC object
#' 
#' Coerce a relatively simple object to an RMINC known object
#' (currently \code{mincList}, \code{mincSingleDim}, \code{mincMultiDim})
#' @param x the object to coerce
#' @return if x is a known minc type return it, if it is a list, attempt
#' to reduce it toa minc object via \link{simplify2minc}, otherwise check
#' if the object has columns, if so reclass it as \code{mincMultiDim} otherwise
#' reclass it as a \code{mincSingleDim}
#' @export
as.minc <- function(x){
  if(is.minc(x)) return(x)
  if(is.list(x)) return(simplify2minc(x))
  
  if(!is.null(ncol(x))){
    class(x) <- c("mincMultiDim", class(x))
    return(x)
  }
  
  class(x) <- c("mincSingleDim", class(x))
  return(x)
} 
    
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
  stopifnot(!is.null(filenames), !is.null(v1))
  enoughAvailableFileDescriptors(length(filenames))
  
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
#' @export
mincGetWorldVoxel <- function(filenames, v1, v2=NULL, v3=NULL) {
  stopifnot(!is.null(filenames), !is.null(v1))
  enoughAvailableFileDescriptors(length(filenames))
  
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
  stopifnot(!is.null(filenames), !is.null(v1), !is.null(v2), !is.null(v3))
  enoughAvailableFileDescriptors(length(filenames))
  
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
  stopifnot(!is.null(filename), !is.null(v1), !is.null(v2), !is.null(v3))
  
  output <- .C("convert_voxel_to_world",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3), PACKAGE="RMINC")$o
  return(output)
}

#' Convert Voxel to World Coordinates
#' 
#' Convert a 3xN matrix of voxel coordinates to world coordinates
#' with respect to a given minc file.
#' 
#' @param filename The minc file indicating a coordinate grid
#' @param voxel_matrix a 3xN matrix of voxel coordinates
#' @return a 3xN matrix of world coordinates
#' @export
mincConvertVoxelMatrix <-
  function(filename, voxel_matrix){
    stopifnot(is.matrix(voxel_matrix), nrow(voxel_matrix) != 3)
    apply(world_matrix, 2, function(row){
      mincConvertWorldToVoxel(filename, row[1], row[2], row[3])
    })
  }

#' World to Voxel
#' 
#' Convert coordinates in the volumes world space (a continuous coordinate space)
#' to a discrete set of voxel coordinates
#' @param filename A path to a minc file
#' @param v1 First world coordinate
#' @param v2 Second world coordinate
#' @param v3 Third world coordinate
#' @param nearest_voxel Logical whether to round the result to the nearest voxel
#' @return a 3-component numeric vector of voxel coordinates
#' @export
mincConvertWorldToVoxel <- function(filename, v1, v2, v3, nearest_voxel = TRUE) {
  stopifnot(!is.null(filename), !is.null(v1), !is.null(v2), !is.null(v3))
  
  output <- .C("convert_world_to_voxel",
               as.character(filename),
               as.double(v1),
               as.double(v2),
               as.double(v3),
               o=double(length=3), PACKAGE="RMINC")$o
  
  if(nearest_voxel)
    output <- round(output)
  
  return(output)
}

#' Convert World to Voxel Coordinates
#' 
#' Convert a 3xN matrix of world coordinates to voxel coordinates
#' with respect to a given minc file.
#' 
#' @param filename The minc file indicating a coordinate grid
#' @param world_matrix a 3xN matrix of world coordinates
#' @param nearest_voxel logical whether to round the results to the nearest voxel
#' @return a 3xN matrix of voxel coordinates
#' @export
mincConvertWorldMatrix <- 
  function(filename, world_matrix, nearest_voxel = TRUE){
    stopifnot(is.matrix(world_matrix), nrow(world_matrix) != 3)
    apply(world_matrix, 2, function(row){
      mincConvertVoxelToWorld(filename, row[1], row[2], row[3], 
                              nearest_voxel = nearest_voxel)
    })
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

#' @method print mincMultiDim
#' @export
print.mincMultiDim <- function(x, ...) {
  cat("Multidimensional MINC volume\n")
  cat("Columns:      ", colnames(x), "\n")
  print(attr(x, "likeVolume"))
}

#' @method print mincLogLikRatio
#' @export 
print.mincLogLikRatio <- function(x, ...) {
  cat("mincLogLikRatio output\n\n")
  mincLmerList <- attr(x, "mincLmerList")
  for (i in 1:length(mincLmerList)) {
    cat("Model", i, ":\n")
    cat("  Formula:  ")
    print(mincLmerList[[i]][[1]]$formula)
  }
  cat("\nMask used:", attr(x, "mask"), "\n")
  cat("Chi-squared Degrees of Freedom:", attr(x, "df"), "\n")
}

#' @method print mincLmer
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
      
#' @method print vertexMultiDim
#' @export
print.vertexMultiDim <- function(x, ...) {
  cat("Multidimensional Vertstats file\n")
  cat("Columns:      ", colnames(x), "\n")
}

#' @method print mincSingleDim     
#' @export 
print.mincSingleDim <- function(x, ...) {
  cat("MINC volume\n")
  print(attr(x, "likeVolume"))
  print(attr(x, "filename"))
}

#' @method print mincQvals
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
  } else {
    file_size <- file.info(as.character(like.filename))$size
    if (is.na(file_size) | length(file_size) == 0 ) {
      stop(c("File ", like.filename, " cannot be found.\n"))
    } 
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
#' @param filename The filename of the MINC2 volume whose dimension sizes are
#' to be returned.
#' @export
minc.dimensions.sizes <- function(filename) {
  stopifnot(!is.null(filename))
  
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

#' @describeIn minc_history retrieve
#' @export
minc.get.history <- 
  function(filename){
    stopifnot(!is.null(filename))
    
    history <- 
      .Call("get_minc_history",
            filename)
    
    unlist(strsplit(history, "\n"))
  }

#' @describeIn minc_history append
#' @export
minc.append.history <-
  function(filename, 
           new_history = NULL){
    
    stopifnot(!is.null(filename))
    
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
  stopifnot(!is.null(filename), !is.null(start), !is.null(count))
  
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
# wilcox.permutation.full <- function(filenames, groupings, mask, n.permute=10) {
#   results <- matrix(nrow=n.permute, ncol=55)
#   mask.volume <- mincGetVolume(mask)
#   for (i in 1:n.permute) {
#     new.order <- sample(groupings)
#     cat("N: ", as.double(new.order), "\n")
#     w <- minc.wilcoxon.test(filenames, new.order, mask)
#     w2 <- w[mask.volume == 1]
#     results[i,] <- tabulate(round(w2),55)
#   }
#   return(results)
# }

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

#' Two tailed version of pt
#' 
#' Generate two tailed probabilities from
#' the cummulative t distribution function with
#' a given number of degrees of freedom.
#' @param q a vector of quantiles
#' @param df degrees of freedom
#' @param log.p whether to return the natural logarithm of p
#' @return A vector of probabilities
#' @seealso pt
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


# create a 2D array of full volumes of all files specified.
minc.get.volumes <- function(filenames) {
  sizes <- minc.dimensions.sizes(filenames[1])
  n.files <- length(filenames)
  output <- matrix(ncol=n.files, nrow=(sizes[1] * sizes[2] * sizes[3]))
  for (i in 1:n.files) {
    output[,i] <- mincGetVolume(filenames[i])
  }
  return(output)
}

#' Create a table of vertex values
#' 
#' Read files containing vertex data into a matrix
#' 
#' @param filenames paths to the vertex data files
#' @return a matrix where each `column` is a matrix of vertex
#' data corresponding to a single file
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
      term <- as.character(formula[[3]])
    }
    else {
      term <- data[[as.character(formula[[3]])]]
    }

    file_covariate <- grepl("\\.mnc$", as.character(term[1]))

    if (file_covariate) {
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
  stopifnot(!is.null(volumeFileName), !is.null(vectorCoord))
  
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

checkCurrentUlimit <- function(){
    
    current_ulimit <- system("ulimit -Sn", intern = TRUE)
    if(current_ulimit == "unlimited") current_ulimit <- Inf
    current_ulimit <- as.numeric(current_ulimit)
    
    return(current_ulimit)
}

tableOpenFiles <- function(){
  lsof_results <- 
    paste("lsof -Ft -p ", Sys.getpid()) %>%
    system(intern = TRUE)
    
  table(lsof_results)  
}

enoughAvailableFileDescriptors <- 
  function(n, error = TRUE){
    available_fds <- checkCurrentUlimit() - sum(tableOpenFiles())
    enough_avail <- (n <= available_fds)
    
    if(error & !enough_avail)
      stop("Not enough available file descriptors to open ",
           n, " files, maybe try setting ulimit -Sn <some-larger-number>")
    
    return(enough_avail)
  }

#' @title Run Testbed
#' @description Run the test bed to ensure all RMINC functions
#' work on your system
#' @param verboseTest
#' Whether or not to verbosely print test output, default is
#' to print simplified results
#' @param dataPath The directory to download and unpack the test data 
#' (unpacks in dataPath/rminctestdata)
#' @param ... additional parameter for \link[testthat]{test_dir}
#' @return invisibly return the test results
#' @export
runRMINCTestbed <- function(..., dataPath = tempdir(), verboseTest = FALSE) {
  
  original_opts <- options()
  old_env_vars <- Sys.getenv(c("TEST_Q_MINC", "NOT_CRAN", "TRAVIS"))
  on.exit({
    options(original_opts)
    do.call(Sys.setenv, as.list(old_env_vars))
  })
  
  setRMINCMaskedValue(0)
  Sys.setenv(TEST_Q_MINC = "yes", NOT_CRAN = "true", TRAVIS = "")

  # Run Tests
  rmincPath <- find.package("RMINC")
  cat("\n\nRunning tests in: ", paste(rmincPath,"/","user_tests/",sep=""), "\n\n\n")
  testReport <- testthat::test_dir(paste(rmincPath,"/","user_tests/",sep=""), ...)
  
  cat("\n*********************************************\n")
  cat("The RMINC test bed finished running all tests\n")
  cat("*********************************************\n\n\n")

  return(invisible(testReport))
}

#' Download Example Data
#' 
#' Download the example data needed to run our examples in your /tmp directory
#' The data can be downloaded manually from 
#' \url{https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz}
#' @param dataPath The directory to download and unpack the test data 
#' (unpacks in dataPath/rminctestdata)
#' @export
getRMINCTestData <- function(dataPath = tempdir()) {

  downloadPath <- file.path(dataPath, "rminctestdata.tar.gz")
  extractedPath <- file.path(dataPath, "rminctestdata/")
  
  if(!file.exists(downloadPath)){
    dir.create(dataPath, showWarnings = FALSE, recursive = TRUE)
    download.file("https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz",
                  destfile = downloadPath,
                  method = "wget") 
  }
  
  untar(downloadPath, exdir = dataPath, compressed = "gzip")
  
  rectifyPaths <-
    function(file){
      readLines(file) %>%
        gsub("/tmp/rminctestdata/", extractedPath, .) %>%
        writeLines(file)
      
      invisible(NULL)
    }
  
  filesToFix <- 
    c("filenames.csv",  
      "minc_summary_test_data.csv",  
      "test_data_set.csv") %>%
    file.path(extractedPath, .)
  
  lapply(filesToFix, rectifyPaths)
  
  invisible(NULL)
}

#' Run function with/without output silenced
#' 
#' used in test bed
#' @param expr an expression to run
#' @param verbose whether to permit output, or capture it
#' @param env the environment in which to run the expression, defaults
#' to the calling environment
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
