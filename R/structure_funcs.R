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
  
