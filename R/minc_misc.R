quiet_mclapply <- function(...){
  sink("/dev/null")
  on.exit(sink())
  
  parallel::mclapply(...)
}

quiet_mcmapply <- function(...){
  sink("/dev/null")
  on.exit(sink())
  
  parallel::mcmapply(...)
}

fix_names <-
  function(x){
    x %>%
      gsub("[^A-Za-z0-9_.]+", "_", .) %>%
      gsub("^_|_$", "", .)
  }

# test to see whether files exist and are readable
minc.isReadable <- function(filenames) {
  rValue <- TRUE
  READ_PERMISSION <- 4
  if (sum(file.access(as.character(filenames), READ_PERMISSION)) != 0
      || is.null(filenames)) {
    
    # ummm, something is not readable, let's say who
    rValue <- FALSE
    for ( ndx in 1:length(filenames) ) {
      if ( file.access(filenames[ndx], READ_PERMISSION) != 0 ) {
        warning(sprintf("File: %s is not accessible.\n"))
      }
    }
  }
  return(rValue)
}


# test to see whether a given file is minc (minc1 or minc2)
isMinc <- function(filename) {
  return(
    file.exists(filename) &
      (isMinc1(filename) | isMinc2(filename)))
}


# test to see whether a given file is minc1
isMinc1 <- function(filename) {
  sysCmd <- paste("file", filename)
  rtnString <- system(sysCmd, intern=TRUE)
  
  return(grepl("NetCDF", rtnString, fixed=TRUE))
}


# test to see whether a given file is minc2
isMinc2 <- function(filename) {
  sysCmd <- paste("file", filename)
  rtnString <- system(sysCmd, intern=TRUE)
  
  return(grepl("Hierachical Data Format", rtnString, fixed=TRUE))
}


# convert minc1 volume to minc2
asMinc2 <- function(filename, output, clobber = FALSE) {
  
  # is it already minc2? Just return the input filename.
  if(isMinc2(filename) ) return(filename)
  
  # if it isn't minc1, tell 'em and run away
  if(!isMinc(filename)) {
    stop(paste("Error: Trying to convert a non-minc file (", filename, ") to minc", sep=""))
  }
  
  # fine. So we now have a minc1 volume that we want to convert to minc2
  #
  # first, get a temporary filename
  cmdOptions <- "-2"
  if(clobber) cmdOptions <- paste(cmdOptions, "-clobber")
  
  # do the conversion
  cat(paste(">> auto-converting", filename, "to minc2 format\n"))
  sysCmd <- paste("mincconvert", cmdOptions, filename, output)
  
  system(sysCmd)
  
  return(invisible(NULL))
}

mincGetDataTypes <- function() {
  # =============================================================================
  # Purpose: return a data.frame containing the minc2 data types
  #
  # Note: Nothing of note, really.
  #
  # =============================================================================
  #
  
  # create a data.frame of data types (an enum in the code)
  # needed to produce human-friendly string
  dataType_numCode <- c(1,3,4,5,6,7,100,101,102,1000,1001,1002,1003,-1)
  dataType_string <- c("8-bit signed integer", "16-bit signed integer", "32-bit signed integer", "32-bit floating point",
                       "64-bit floating point", "ASCII string", "8-bit unsigned integer", "16-bit unsigned integer", "32-bit unsigned integer",
                       "16-bit signed integer complex", "32-bit signed integer complex", "32-bit floating point complex",
                       "64-bit floating point complex", "non_uniform record")
  dataType_code <- c("MI_TYPE_BYTE", "MI_TYPE_SHORT", "MI_TYPE_INT", "MI_TYPE_FLOAT", "MI_TYPE_DOUBLE",
                     "MI_TYPE_STRING", "MI_TYPE_UBYTE", "MI_TYPE_USHORT", "MI_TYPE_UINT", "MI_TYPE_SCOMPLEX",
                     "MI_TYPE_ICOMPLEX", "MI_TYPE_FCOMPLEX", "MI_TYPE_DCOMPLEX", "MI_TYPE_UNKNOWN")
  dataTypes.df <- data.frame(code=dataType_code, numCode=dataType_numCode, string=dataType_string, stringsAsFactors=FALSE)
  
  # return the valid data types
  return(dataTypes.df)
}



mincGetDataClasses <- function() {
  # =============================================================================
  # Purpose: return a data.frame containing the minc2 data classes
  #
  # Note: Nothing of note, really.
  #
  # =============================================================================
  #
  
  # create a data.frame of data classes (an enum in the code)
  # needed to produce human-friendly string
  dataClass_numCode <- c(0, 1, 2, 3, 4, 5)
  dataClass_string <- c("REAL", "INTEGER", "LABEL", "COMPLEX", "UNIFORM_RECORD", "NON_UNIFORM_RECORD")
  dataClasses.df <- data.frame(numCode=dataClass_numCode, string=dataClass_string, stringsAsFactors=FALSE)
  
  # return the valid data classes
  return(dataClasses.df)
}


rminc.readLinearXfmFile <- function(xfmFilename) {
  program <- "xfm2param"
  progOptions <- "-version"
  test_string <- "mni_autoreg"
  status <- rminc.checkForExternalProgram(program, test_string, progOptions)
  if ( !status ) { stop("Program xfm2param of package mni_autoreg cannot be found on PATH") }
  # OK, so we have xfm2pram -- now let's do the read
  # ... first have xfm2param place tabular output in a temp file
  tmpfile <- tempfile("rminc.readLinearXfmFile")
  cmd <- paste("xfm2param", xfmFilename, "> ", tmpfile)
  # ... now read the nicely formatted tabular file
  system(cmd, intern=TRUE, wait=TRUE)
  xfm.df <- read.table(tmpfile, skip=1, stringsAsFactors=FALSE)
  
  # make first column into row names
  rowNames <- xfm.df[,1]
  rowNames <- gsub("^-", "", rowNames)
  row.names(xfm.df) <- rowNames
  
  # change col names and then remove col 1
  names(xfm.df) <- c("dummy", "x", "y", "z")
  xfm.df <- subset(xfm.df,select=-dummy) 
  
  # return the xfm data.frame
  return(xfm.df)
}

rminc.checkForExternalProgram <- function(program, test_string, prog_options="") {
  cmd <- paste(program, prog_options)
  cmdOut <- system(cmd, intern=TRUE, wait=TRUE)
  
  # collapse all output into a single line for easy grepping
  cmdOutLong <- paste(cmdOut, collapse="")
  
  # look for test string in output
  if ( !grepl(test_string, cmdOutLong, fixed=TRUE) ) {
    # test string not found??
    cat(sprintf("Attempt to execute program \"%s\" within shell failed\n", program))
    cat(sprintf("Shell responded with: \n%s\n", cmdOut))
    warning("\nCheck your path ...")
    return(FALSE)
  } 
  
  # return TRUE if we made this far
  return(TRUE)
}



#' Explode a Label Volume into its Components
#' 
#' Given a label volume, this function splits the volume by label and then
#' returns a list() containing a mask volume for each of the labels.
#' 
#' @param label_vol A string containing the fully-qualified path to the input
#' label volume.
#' @param labels Options vector of label names
#' @param civetLabels A logical variable indicating whether the label volume is
#' using the Civet convention with regards to naming tissue type (e.g.,
#' 0=background, 1=csf, etc). If TRUE, the returned list components are named
#' using Civet tissue types (bg, csf, gm. wm), else components are simply
#' labelled by label number e.g. (``label_0'', ``label_2'', etc.).
#' @return A list is returned with each list item holding a mask volume
#' reflecting a particular label.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
volume.explodeLabelVolume <- function(label_vol, labels=NULL, civetLabels=TRUE) {
  
  # list the tissue types in Civet order, to be used in naming the
  # components of the returned list
  tissueTypes <- c("bg", "csf", "gm", "wm")
  
  # read the label volume
  # label_vol <- mincIO.readVolume(labelVolname)
  
  # here we have a choice: We either specify which labels we want exploded out,
  # or, we don't, and we get all of them
  #
  if ( is.null(labels)) {
    # not specified, so get the set of unique labels
    labels <- unique(label_vol)
  } else {
    # specified, ... make sure we've been passed a numeric vector
    if ( !is.vector(labels, mode="numeric")) {stop("Label vector must be a numeric vector")}
    
  }
  
  
  # loop thru all labels, creating a mask for each
  return_list <- list()
  for ( label in labels ) {
    #
    # determine the list component name
    label_name <- paste("label", label, sep="_")
    if ( civetLabels ) { label_name <- tissueTypes[label +1]}
    
    # compute the mask
    #cat(sprintf("processing label %s\n", label_name))
    mask_vol <- ifelse(label_vol == label, 1, 0)
    
    # adjust volume attributes to type of "mask"
    # mincIO.setProperty(mask_vol, "volumeType", "mask")
    # mincIO.setProperty(mask_vol, "colorMap", "gray")
    # mincIO.setProperty(mask_vol, "volumeIntensityRange", c(0,1))
    
    # save it to the return list
    return_list[[label_name]] <- mask_vol
  }
  
  return(return_list)
}


#' Combine Multiple Mask Volumes into a Single Mask Volume 
#' 
#' Given a list containing more than one label volume, combine those 
#' volumes, creating an aggregate mask volume. 
#' 
#' @param vol_list A list containing more than one mask volume.  Note that all 
#' volumes must reflect the same sampling. 
#' @return A single aggregate MincVolumeIO volume is returned. 
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
volume.combineMaskVolumes <- function(vol_list) {
  
  # first, get the number of volumes to combine
  nVolumes <- length(vol_list)
  
  # loop thru all labels, adding 'em up
  vol_sum <- vol_list[[1]]
  for ( vol in 2:nVolumes ) {
    vol_sum <- vol_sum + vol_list[[vol]]
  }
  
  # set voxels included in multiple masks to '1'
  vol_sum <- ifelse(round(vol_sum) < 1, 0, 1)
  #
  # make sure it's still a volume, then return
  #vol_sum <- mincIO.asVolume(vol_sum, vol_list[[1]])
  return(vol_sum)
}

# a replica of the qvalue function from the "qvalue" package, but with
# certain parts modified (see comments in code) for speed.
# Should produce identical results.

fast.qvalue <- function (p, lambda = seq(0, 0.9, 0.05),
                         pi0.method = "smoother", 
                         fdr.level = NULL, robust = FALSE,
                         gui = FALSE, smooth.df = 3, 
                         smooth.log.pi0 = FALSE) 
{
  if (gui & !interactive()) 
    gui = FALSE
  if (min(p) < 0 || max(p) > 1) {
    if (gui) 
      eval(expression(postMsg(paste("ERROR: p-values not in valid range.", 
                                    "\n"))), parent.frame())
    else print("ERROR: p-values not in valid range.")
    return(0)
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    if (gui) 
      eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.", 
                                    "\n"))), parent.frame())
    else print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
    return(0)
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 
                             1)) {
    if (gui) 
      eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                                    "\n"))), parent.frame())
    else print("ERROR: Lambda must be within [0, 1).")
    return(0)
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      if (gui) 
        eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", 
                                      "\n"))), parent.frame())
      else print("ERROR: Lambda must be within [0, 1).")
      return(0)
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
  }
  else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) 
        pi0 <- log(pi0)
      spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
      pi0 <- predict(spi0, x = max(lambda))$y
      if (smooth.log.pi0) 
        pi0 <- exp(pi0)
      pi0 <- min(pi0, 1)
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - 
                                                     lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
    }
    else {
      print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
      return(0)
    }
  }
  if (pi0 <= 0) {
    if (gui) 
      eval(expression(postMsg(paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.", 
                                    "\n"))), parent.frame())
    else print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
    return(0)
  }
  if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
                              1)) {
    if (gui) 
      eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", 
                                    "\n"))), parent.frame())
    else print("ERROR: 'fdr.level' must be within (0, 1].")
    return(0)
  }
  u <- order(p)
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  if (robust) {
    qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
  }
  #qvalue[u[m]] <- min(qvalue[u[m]], 1)
  
  qvalue <- .C("qvalue_min", as.double(qvalue), as.integer(u),
               as.integer(m),
               o=double(length=m), PACKAGE="RMINC")$o
  
  #for (i in (m - 1):1) {
  #    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  #}
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
                   pvalues = p, fdr.level = fdr.level, significant = (qvalue <= 
                                                                        fdr.level), lambda = lambda)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
                   pvalues = p, lambda = lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}



