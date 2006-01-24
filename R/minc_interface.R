minc.get.voxel.from.files <- function(filenames, v1, v2, v3) {
  num.files <- length(filenames)
  output <- .C("get_voxel_from_files",
               as.character(filenames),
               as.integer(num.files),
               as.integer(v1),
               as.integer(v2),
               as.integer(v3),
               o=double(length=num.files))$o
  return(output)
}


minc.get.volume <- function(filename) {
  sizes <- minc.dimensions.sizes(filename)
  start <- c(0,0,0)
  total.size <- sizes[1] * sizes[2] * sizes[3]
  output <- .C("get_hyperslab",
               as.character(filename),
               as.integer(start),
               as.integer(sizes),
               hs=double(total.size))$hs
  return(output)
}

minc.write.volume <- function(output.filename, like.filename, buffer) {
  sizes <- minc.dimensions.sizes(like.filename)
  start <- c(0,0,0)
  if ((sizes[1] * sizes[2] * sizes[3]) != length(buffer)) {
    error("Size of like-file not the same as size of buffer")
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
  

# given an array of filenames, get the data from the corresponding slices
minc.get.slices <- function(filenames, begin.slice=NA, end.slice=NA) {
  # blithely assume all volumes have the same dimensions
  #cat(filenames[1]
  sizes <- .C("get_volume_sizes",
              as.character(filenames[1]),
              sizes = integer(3))$sizes
  num.slices <- end.slice - begin.slice
  total.size <- sizes[2] * sizes[3] * num.slices
  num.files <- length(filenames)
  start <- c(begin.slice, 0, 0)
  count <- c((end.slice - begin.slice), sizes[2], sizes[3])
  output <- matrix(ncol=num.files, nrow=total.size)

  # get the hyperslabs
  for (i in 1:num.files) {
    output[, i] <- .C("get_hyperslab",
                       as.character(filenames[i]),
                       as.integer(start),
                       as.integer(count), hs=double(total.size))$hs
  }
  return(output)
}

minc.dimensions.sizes <- function(filename) {
  sizes <- .C("get_volume_sizes",
              as.character(filename),
              sizes = integer(3))$sizes
  return(sizes)
}

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

minc.get.hyperslab2 <- function(filename, start, count, buffer=NA) {
#  total.size <- (count[1] - start[1]) * (count[2] - start[2]) *
#                (count[3] - start[3])
  total.size <- (count[1] * count[2] * count[3])
  if (is.na(buffer)) {
    buffer <- double(total.size)
  }
  cat("Total size: "); cat(total.size); cat("\n")
  .Call("get_hyperslab2",
               as.character(filename),
               as.integer(start),
               as.integer(count), hs=buffer, DUP=FALSE)
  #return(output)
}

wilcox.permutation <- function(filenames, groupings, mask, n.permute=10) {
  results <- data.frame(greater50=vector(length=n.permute),
                        less5=vector(length=n.permute),
                        equal55=vector(length=n.permute),
                        equal0=vector(length=n.permute),
                        percent5=vector(length=n.permute),
                        percent95=vector(length=n.permute))
  mask.volume <- minc.get.volume(mask)
  for (i in 1:n.permute) {
    new.order <- sample(groupings)
    cat("N: ", as.double(new.order), "\n")
    w <- minc.wilcoxon.test(filenames, new.order, mask)
    w2 <- w[mask.volume == 1]
    ws <- sort(w2)
    #minc.write.volume(paste("p", i, ".mnc", sep=""), mask, w)
    results$greater50[i] <- sum(w2 > 50)
    results$less5[i] <- sum(w2 < 5)
    results$equal55[i] <- sum(w2 == 55)
    results$equal0[i] <- sum(w2 == 0)
    results$percent5[i] <- ws[round(length(ws) * 0.05)]
    results$percent95[i] <- ws[round(length(ws) * 0.95)]
  }
  return(results)
}

minc.wilcoxon.test <- function(filenames, groupings, mask=NULL) {
  voxel.wilcoxon.test <- function(x) {
    .Call("wilcoxon_rank_test", as.double(x),
          as.double(tmp.groupings),
          as.double(sum(tmp.groupings==0)),
          as.double(sum(tmp.groupings==1)))
  }

  groupings <- as.double(groupings)
  assign("tmp.groupings", groupings, env=.GlobalEnv)
  cat("TMP: ", as.double(tmp.groupings), "\n")
  assign("voxel.wilcoxon.test", voxel.wilcoxon.test, env=.GlobalEnv)
  w <- minc.apply(filenames, quote(voxel.wilcoxon.test(x)), mask)
  return(w)
}

minc.read.volumes <- function(filenames) {
  sizes <- minc.dimensions.sizes(filenames[1])
  n.files <- length(filenames)
  output <- matrix(ncol=n.files, nrow=(sizes[1] * sizes[2] * sizes[3]))
  for (i in 1:n.files) {
    output[,i] <- minc.get.volume(filenames[i])
  }
  return(output)
}

# efficient way of applying an R function to every voxel
minc.apply <- function(filenames, function.string, mask=NULL) {
  .Call("minc2_apply", as.character(filenames),
        function.string,
        as.double(! is.null(mask)),
        as.character(mask),
        parent.env(environment()))
}

# perform specified function across all voxels of all volumes
minc.slice.loop <- function(filenames, n.slices, output.dims, func, ...) {
  sizes <- minc.dimensions.sizes(filenames[1])
  slice.dim.size <- sizes[1]

  total.size <- sizes[1] * sizes[2] * sizes[3]
  n.subjects <- length(filenames)

  sequence <- seq(0, slice.dim.size, n.slices)

  # make sure that all slices are included
  if (sequence[length(sequence)] != slice.dim.size) {
    sequence <- append(sequence, slice.dim.size)
  }
  cat("Sequence: "); cat(sequence); cat("\n");

  output <- vector(length=total.size)
  # do the loop
  prev.index <- 1
  buffer <- matrix(0, nrow=((sequence[2] - sequence[1]) *sizes[2] * sizes[3]), ncol=n.subjects)

  for (i in 1:(length(sequence)-1)) {
    start <- c(sequence[i], 0, 0)
    count <- c((sequence[i+1] - sequence[i]), sizes[2], sizes[3])

    #start <- c(0,0, sequence[i])
    #count <- c(sizes[1], sizes[2], sequence[i+1] - sequence[i])
    
    current.index <- sequence[i+1] * sizes[2] * sizes[3]
    #current.size <- current.index - prev.index
    current.size <- count[1] * count[2] * count[3]
    cat("Counter: "); cat(i); cat("\n")
    cat("sequence: "); cat(sequence[i]); cat("\n")
    cat("next sequence: "); cat(sequence[i+1]); cat("\n")
    cat("start: "); cat(start); cat("\n")
    cat("count: "); cat(count); cat("\n")
    for (j in 1:n.subjects) {
      #buffer[,j] <- minc.get.hyperslab2(filenames[j], start, count, buffer[,j])
    }
    current.index <- current.index
    cat("Prev: "); cat(prev.index); cat(" Curr: ");
    cat(current.index); cat("\n")
    for (k in prev.index:current.index) {
      output[k] <- wt(buffer[k,])
    }
    #output[prev.index:current.index] <- apply(buffer, 1, func, ...)
    #output[prev.index:current.index] <- buffer[,1]
    prev.index <- prev.index + current.size
    #sequence[i+1] <- sequence[i+1] + 1
  }
  return(output)
}
