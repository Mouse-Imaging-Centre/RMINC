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
  output <- system(paste(command, collapse=" "), intern = TRUE)
  cat(output, sep = "\n")
  
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
  # make a data frame out of the tags
  tags <- as.data.frame(tags)
  return(tags)
}


#' read a tag file
#'
#' @param filename string containing filename to read 
#'
#' @return tags as a ntags by the number of fields in the input tag matrix
#' @export
mincGetTagFile <- function(filename) {
  tags <- readLines(filename)
  # tag points begin after Points = line
  beginIndex <- grep("Points", tags) + 1
  tags <- tags[beginIndex:length(tags)]
  # get rid of trailing semicolon
  tags[length(tags)] <- sub(";", "", tags[length(tags)]) 
  # get rid of starting whitespace
  tags <- sub("^ *", "", tags)
  coord_table <- read.table(text = tags)
  coord_table %>% as.matrix() %>% unname
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

#' label peaks with the name of the atlas structure they are in
#'
#' @param peaks the output of \code{\link{mincFindPeaks}}
#' @param atlas the atlas volume, either as a \code{\link{mincArray}} or as a filename
#' @param defs the atlas definitions, same as used for \code{\link{anatGetAll}}
#'
#' @return the peaks data frame, with an extra column containing the label
#' @export
#'
#' @examples
#' \dontrun{
#' peaks <- mincFindPeaks(-log10(qvs), "Neonatal:time.to.sac", 
#'                        "pos", posThreshold=1.3, minDistance=1)
#' peaks <- mincLabelPeaks(peaks, 
#'                          atlasVol, 
#'                          defs="Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv")
#' }
mincLabelPeaks <- function(peaks, atlas, defs=getOption("RMINC_LABEL_DEFINITIONS")) {
  
  # if atlas is a filename, then read it in as a mincArray
  if (is.character(atlas)) {
    atlas <- mincArray(mincGetVolume(atlas))
  }
  else if (length(dim(atlas)) == 3) {
    # do nothing
  }
  else {
    stop("Error: atlas has to be either a mincArray or a filename")
  }
  
  # get the values from the atlas at every index of the peaks
  peaks$label <- atlas[as.matrix(peaks[,1:3])]
  # defs is not set to NULL, then translate numbers to label names
  if (!is.null(defs)) {
    defs <- read.csv(defs)
    # melt into long format
    mdefs <- defs[,1:3] %>% gather_("variable", "value", c("right.label", "left.label"))
    mdefs$variable <- as.character(mdefs$variable)
    mdefs$Structure <- as.character(mdefs$Structure)
    # get rid of .label
    mdefs$variable <- sub('.label', '', mdefs$variable, fixed=T)
    # change the side to empty if label is bilateral
    dups <- mdefs[duplicated(mdefs$value),"value"]
    for (i in 1:length(dups)) {
      mdefs[mdefs$value == dups[i],"variable"] <- ""
    }
    # merge side and structure name
    mdefs$Structure <- paste(mdefs[,2], mdefs[,1])
    # add an unlabelled structure for peaks outside of the label volume
    mdefs <- rbind(mdefs, c("unlabelled", "both", 0))
    # now set structure names
    for (i in 1:nrow(peaks)) {
      peaks$label[i] <- mdefs$Structure[mdefs$value == peaks$label[i]]
    }
  }
  return(peaks)
}
