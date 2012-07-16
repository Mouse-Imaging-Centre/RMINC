# compute a linear model over every structure


anatGetFile <- function(filename, atlas, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv", dropLabels=FALSE, side="both" ) {
  out <- NULL
  if (method == "jacobians") {
    system(paste("label_volumes_from_jacobians", atlas, filename, "> tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "labels") {
    # filename here should be a set of labels unique to this brain
    system(paste("volumes_from_labels_only.py", filename, "tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "means") {
    system(paste("average_values_across_segmentation.py", filename, atlas, "tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "text") {
    # values are already extracted and stored in a text file
    out <- read.table(filename, header=FALSE)
  }
  #cat("FILENAME:", filename, "\n")
  if (dropLabels == TRUE) {
    labels <- read.csv(defs)
    usedlabels <- 0
    if (side == "right") {
      usedlabels <- labels$right.label
    }
    else if (side == "left") {
      usedlabels <- labels$left.label
    }
    else if (side == "both") {
      usedlabels <- c(labels$left.label, labels$right.label)
    }
    out <- out[out$V1 %in% usedlabels,]
  }
  return(out)
}

anatRenameRows <- function(anat, defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv") {
  defs <- read.csv(defs)
  rn <- rownames(anat)
  on <- as.character(rn)
  for (i in 1:length(rn)) {
    # if structure exists in both hemisphere with same number
    if (rn[i] %in% defs$right.label & rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste(defs$Structure[index])
      #cat("Match", rn[i], "to both:", on[i], "\n")
    }
    # left side only
    else if (rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$left.label)
      on[i] <- paste("left", defs$Structure[index])
      #cat("Match", rn[i], "to left:", on[i], "\n")
    }
    else if (rn[i] %in% defs$right.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste("right", defs$Structure[index])
      #cat("Match", rn[i], "to right:", on[i], "\n")
    }
    else {
      #cat("Something weird", rn[i], "\n")
    }
  }
  rownames(anat) <- on
  # make it an anat class
  class(anat) <- c("anatUnilateral", "anatMatrix", "matrix")
  # make sure we can map label numbers for later use
  attr(anat, "anatIDs") <- rn
  return(anat)
}

# both unilateral and bilateral matrices can be printed the same way
print.anatMatrix <- function(x, ...) {
  print.table(x)
}

anatGetAll <- function(filenames, atlas, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv", dropLabels=FALSE, side="both") {
  vol <- anatGetFile(filenames[1], atlas, method, defs, dropLabels, side)
  output <- matrix(nrow=nrow(vol), ncol=length(filenames))
  rownames(output) <- vol[,1]

  output[,1] <- vol[,2]
  for (i in 2:length(filenames)) {
    output[,i] <- anatGetFile(filenames[i], atlas, method, defs, dropLabels, side)[,2]
  }
  attr(output, "atlas") <- atlas
  if (! is.null(defs)) {
    output <- anatRenameRows(output, defs)
  }
  return(t(output))
}

anatCombineStructures <- function(vols, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv") {
  labels <- read.csv(defs)
  combined.labels <- matrix(nrow=nrow(vols), ncol=nrow(labels))
  labelNumbers <- attr(vols, "anatIDs")
  for (i in 1:nrow(labels)) {
    if (labels$right.label[i] == labels$left.label[i]) {
      combined.labels[,i] <- vols[, labelNumbers == as.character(labels$right.label[i])]
    }
    else {
      if (method == "jacobians" || method == "labels") {
        combined.labels[,i] <-
          vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]
      }
      else if (method == "means"){
        combined.labels[,i] <-
          (vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]) /2
      }
        
    }
  }
  colnames(combined.labels) <- labels$Structure
  class(combined.labels) <- c("anatCombined", "anatMatrix", "matrix")
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
  result <- .Call("vertex_lm_loop", anatmatrix, mmatrix, PACKAGE="RMINC")
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
  attr(result, "stat-type") <- c("F", rep("t", ncol(result)-1))
  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  dflist <- vector("list", ncol(result))
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  
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
