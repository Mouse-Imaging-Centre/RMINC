
#' Finds peak values in a MINC file or volume
#' 
#' Inputs a MINC file or RMINC output and identifies peak voxels (given some constraints) within
#'
#' @param inputStats Either an RMINC object or a filename from which peaks are to be determined.
#' @param column If an RMINC object then the column name or number
#' @param direction "both", "pos" for positive peaks only, or "neg" for negative peaks only
#' @param minDistance minimum distance, in mm, between peaks
#' @param threshold threshold above which to consider peaks
#' @param posThreshold positive threshold (overrides threshold if both are given)
#' @param negThreshold negative threshold (overrides threshold if both are given)
#' 
#' @details Provides an interface to find_peaks, a command line program that comes with the MINC library.
#'
#' @return The set of peaks as a matrix. The matrix is of dimensions ntags x 7; the first
#' three columns correspond to the first three dimensions in the same coordinates as used by 
#' mincArray; the next three columns are the x, y, and z coordinates in world space, and the
#' 7th column is the value at that location.
#' @export
mincFindPeaks <- function(inputStats, column=1, direction="both", minDistance=NA,
                          threshold=NA, posThreshold=NA, negThreshold=NA) {
  
  # test if input is a file or the output of mincLm and friends
  if (is.character(inputStats)) {
    filename = inputStats
  }
  else {
    filename = tempfile(fileext = ".mnc")
    mincWriteVolume(inputStats, filename, column)
  }
  
  # output tags as a temporary file
  tagfile = tempfile(fileext = ".tag")
  
  # hold the command in a vector
  command <- c("find_peaks", filename, tagfile)
  
  # set directions
  if( !(direction %in% c("both", "pos", "neg"))) {
    stop("direction argument must be one of both, pos, or neg")
  }
  
  if (direction == "pos") {
    command <- append(command, "-pos_only")
    
  }
  else if (direction == "neg") {
    command <- append(command, "-neg_only")
  }
  
  # set thresholds
  if (!all(is.na(c(threshold, posThreshold)))) {
    if (is.na(posThreshold)) {
      posThreshold <- threshold
    }
    command <- append(command, c("-pos_threshold", posThreshold))
  }
  
  if (!all(is.na(c(threshold, negThreshold)))) {
    if (is.na(negThreshold)) {
      negThreshold <- -threshold
    }
    command <- append(command, c("-neg_threshold", negThreshold))
  } 
  
  if (!is.na(minDistance)) {
    command <- append(command, c("-min_distance", minDistance))
  }
  
  # call find_peaks on command line
  system(paste(command, collapse=" "))
  
  # read in the resulting tags
  tags <- mincGetTagFile(tagfile)
  # and convert to mincArray coordinates
  tagsMA <- mincConvertTagToMincArrayCoordinates(tags, filename)
  # merge the two together
  tags <- cbind(tagsMA[,1:3], tags[,c(1:3,7)])
  colnames(tags) <- c("d1", "d2", "d3", "x", "y", "z", "value")
  
  # delete temporary files
  file.remove(tagfile)
  if (!is.character(inputStats)) {
    file.remove(filename)
  }
  return(tags)
}


#' read a tag file
#'
#' @param filename string containing filename to read 
#'
#' @return tags as a ntags by 7 matrix
#' @export
mincGetTagFile <- function(filename) {
  tags <- scan(filename, what="character")
  # tag points begin after Points = line
  beginIndex <- grep("Points", tags) + 2
  endIndex <- length(tags)-1 # get rid of training ;
  return(matrix(as.numeric(tags[beginIndex:endIndex]), ncol=7, byrow=T))
}

#' convert tags from world coordinates to mincArray coordinates
#'
#' @param tags from mincGetTagFile
#' @param filename path to a MINC file to obtain coordinate info from
#'
#' @return tags as a ntags by 4 matrix
#' @export
mincConvertTagToMincArrayCoordinates <- function(tags, filename) {
  tags <- tags[,c(1:3,7)] # get rid of repeated coordinates
  out <- matrix(ncol=ncol(tags), nrow=nrow(tags))
  for (i in 1:nrow(tags)) {
    out[i,3:1] <- mincConvertWorldToVoxel(filename, tags[i,1], tags[i,2], tags[i,3]) + 1
  }
  out[,4] <- tags[,4]
  return(out)
}

