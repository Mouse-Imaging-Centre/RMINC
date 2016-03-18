
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







