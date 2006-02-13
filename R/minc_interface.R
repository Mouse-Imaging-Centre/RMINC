# get the real value of one voxel from all files.
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

# return a volume as a 1D array.
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

# write a 1D buffer into a new file, getting volume info from an
# existing file.
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

# get the dimension sizes of a particular file.
minc.dimensions.sizes <- function(filename) {
  sizes <- .C("get_volume_sizes",
              as.character(filename),
              sizes = integer(3))$sizes
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
               as.integer(count), hs=buffer)$hs
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

# run either a t-test or wilcoxon test at every voxel
minc.model <- function(filenames, groupings, method="t-test",
                       mask=NULL) {

  if (method == "t-test" || method == "wilcoxon"
      || method == "correlation" || method == "lm") {
    # do nothing
  }
  else {
    stop("Method must be one of t-test, wilcoxon, correlation or lm")
  }

  if (method == "lm") {
    .Call("minc2_model",
          as.character(filenames),
          as.matrix(groupings),
          as.double(! is.null(mask)),
          as.character(mask),
          as.character(method))
  }
  else {
    groupings <- as.double(groupings)
    .Call("minc2_model",
          as.character(filenames),
          as.double(groupings),
          as.double(! is.null(mask)),
          as.character(mask),
          as.character(method))
  }
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

# efficient way of applying an R function to every voxel
minc.apply <- function(filenames, function.string, mask=NULL) {
  .Call("minc2_apply", as.character(filenames),
        function.string,
        as.double(! is.null(mask)),
        as.character(mask),
        parent.env(environment()))
}

