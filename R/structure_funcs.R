# compute a linear model over every structure


anatGetFile <- function(filename, atlas, method="volume") {
  if (method == "volume") {
    system(paste("label_volumes_from_jacobians", atlas, filename, "> tmp.txt", sep=" "))
  }
  else if (method == "means") {
    system(paste("average_values_across_segmentation.py", filename, atlas, "tmp.txt", sep=" "))
  }
  return(read.csv("tmp.txt", header=F))
}

anatGetAll <- function(filenames, atlas, method="volume") {
  vol <- anatGetFile(filenames[1], atlas, method)
  rownames(vol) <- vol[,1]
  vol[,1] <- vol[,2]
  for (i in 2:length(filenames)) {
    vol[,i] <- getvol(filenames[i], atlas, method)[,2]
  }
  attr(vol, "atlas") <- atlas
  return(vol)
}

anatCombineStructures <- function(vols, method="volume", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv") {
  labels <- read.csv(defs)
  combined.labels <- matrix(nrow=ncol(vols), ncol=nrow(labels))
  for (i in 1:nrow(labels)) {
    if (labels$right.label[i] == labels$left.label[i]) {
      combined.labels[,i] <- t(vols[as.character(labels$right.label[i]),])
    }
    else {
      if (method == "volume") {
        combined.labels[,i] <- t(vols[as.character(labels$right.label[i]),] +
                                 vols[as.character(labels$left.label[i]),])
      }
      else {
        combined.labels[,i] <- t(vols[as.character(labels$right.label[i]),] +
                                 vols[as.character(labels$left.label[i]),]) /2
      }
        
    }
  }
  colnames(combined.labels) <- labels$Structure
  class(combined.labels) <- c("anatMatrix", "matrix")
  attr(combined.labels, "atlas") <- attr(vols, "atlas")
  attr(combined.labels, "definitions") <- defs
  return(combined.labels)
}

anatApply <- function(vols, grouping, method=mean) {
  ngroups <- length(levels(grouping))
  output <- matrix(nrow=ncol(vols), ncol=ngroups)

  for (i in 1:ncol(vols)) {
    output[i,] <- tapply(vols[,i], grouping, method)
  }
  colnames(output) <- levels(grouping)
  rownames(output) <- colnames(vols)
  return(output)
}
  
anatLm <- function(formula, data, anat, subset=NULL) {
  # the same code to extract the formula as in mincLm ...
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  # ... except this time we add a row to keep track of what subsetting does
  mf$rowcount <- seq(1, nrow(data))
  m <- match(c("formula", "data", "subset", "rowcount"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  mmatrix <- model.matrix(formula, mf)
  # get the data from the anatomy matrix using the same subset as specified
  # in the formulat
  anatmatrix <- t(anat[mf[,"(rowcount)"],])
  # same stats as for vertex tables
  result <- .Call("vertex_lm_loop", anatmatrix, mmatrix)
  rownames(result) <- colnames(anat)
  # get the first voxel in order to get the dimension names
  v.firstVoxel <- anatmatrix[1,]
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  colnames(result) <- c("F-statistic", rows)
  class(result) <- c("anatModel", "matrix")
  attr(result, "atlas") <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model") <- as.matrix(mmatrix)
  return(result)
    
}

anatCreateVolume <- function(anat, filename, column=1) {
  labels <- read.csv(attr(anat, "definitions"))
  volume <- mincGetVolume(attr(anat, "atlas"))
  newvolume <- volume
  for (i in 1:nrow(labels)) {
    newvolume[volume < labels$right.label[i] + 0.5 &
              volume > labels$right.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
    newvolume[volume < labels$left.label[i] + 0.5 &
              volume > labels$left.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
  }
  mincWriteVolume(newvolume, filename, attr(anat, "atlas"))
}

anatFDR <- function(buffer, method="FDR") {
  vertexFDR(buffer, method)
}
