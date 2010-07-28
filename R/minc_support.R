
# test to see whether files exist and are readable
rminc.isReadable <- function(filenames) {
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
rminc.isMinc <- function(filename) {
	rValue <- FALSE
	if ( !file.exists(filename) ) { return(rValue) }
	if ( rminc.isMinc1(filename) ) rValue <- TRUE
	if ( rminc.isMinc2(filename) ) rValue <- TRUE
	#
	return(rValue)
}


# test to see whether a given file is minc1
rminc.isMinc1 <- function(filename) {
	rValue <- FALSE
	sysCmd <- paste("file", filename)
#	print(sysCmd)
	rtnString <- system(sysCmd, intern=TRUE)
#	print(rtnString)
	if ( grepl("NetCDF", rtnString, fixed=TRUE) ) rValue <- TRUE
	return(rValue)
}


# test to see whether a given file is minc2
rminc.isMinc2 <- function(filename) {
	rValue <- FALSE
	sysCmd <- paste("file", filename)
#	print(sysCmd)
	rtnString <- system(sysCmd, intern=TRUE)
#	print(rtnString)
	if ( grepl("Hierarchical Data Format", rtnString, fixed=TRUE) ) rValue <- TRUE
	return(rValue)
}


# convert minc1 volume to minc2
rminc.asMinc2 <- function(filename, keepName=TRUE) {
	
	# is it already minc2? Just return the input filename.
	if ( rminc.isMinc2(filename) ) return(filename)
	
	# if it isn't minc1, tell 'em and run away
	if ( !rminc.isMinc(filename) ) {
		stop(paste("Error: Trying to convert a non-minc file (", filename, ") to minc", sep=""))
	}
	
	# fine. So we now have a minc1 volume that we want to convert to minc2
	#
	# first, get a temporary filename
	cmdOptions <- ""
	if ( keepName ) {
		# we want to use the input filename, but put the file in tmpdir
		# ... allow for overwrite of file in tmpdir
		cmdOptions <- "-clobber"
		tmpFile <- basename(filename)
		tmpFile <- file.path(tempdir(), tmpFile)
		}
	else {
		tmpFile <- tempfile( pattern="R_mincIO_mincconvert_")
	}
	
	# do the conversion
	cat(paste(">> auto-converting", filename, "to minc2 format\n"))
	sysCmd <- paste("mincconvert", cmdOptions, "-2",  filename, tmpFile)
#	print(sysCmd)
	system(sysCmd, wait=TRUE)
	#
	return(tmpFile)
}


rminc.convertVoxelToWorld <- function(filename, voxCoords) {
	#
	
	if ( R_DEBUG_mincIO ) cat(sprintf(">>rminc.convertVoxelToWorld\n"))

	# the C routines want the coordinates 0-relative, AND in volume order
	# ... so convert first
	voxCoords <- rev(voxCoords) -1
	
	# dunno why, but the voxel coordinates are passed as doubles
	output <- .Call("convert_voxel_to_world_mincIO",
               as.character(filename),
               as.double(voxCoords), PACKAGE="RMINC")

	# return a vector of 3 doubles
	if ( R_DEBUG_mincIO ) cat(sprintf("<<rminc.convertVoxelToWorld\n"))
	return(output)
}



rminc.convertWorldToVoxel <- function(filename, worldCoords) {
	#
	# dunno why, but the voxel coordinates are passed as doubles
	if ( R_DEBUG_mincIO ) cat(sprintf(">>rminc.convertWorldToVoxel\n"))

	output <- .Call("convert_world_to_voxel_mincIO",
               as.character(filename),
               as.double(worldCoords), PACKAGE="RMINC")
	if ( R_DEBUG_mincIO ) {
	  cat(sprintf("voxel coordinates returned by .Call ...\n"))
	  print(output)
	}

	# return a vector of 3 doubles 
	# ... 0-relative and (t),z,y,x order from C, so let's adjust it for R
	if ( length(output) == 3 ) {
		output <- rev(output) +1
		
	} else {
		# 4D volume, so ignore the (first) time dimension
		output <- rev(output[2:4]) +1
	}
	
	# send it back
	if ( R_DEBUG_mincIO ) cat(sprintf("<<rminc.convertWorldToVoxel\n"))
	return(output)
}



rminc.getDataTypes <- function() {
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



rminc.getDataClasses <- function() {
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



rminc.checkForExternalProgram <- function(program, test_string, prog_options="") {
	# =============================================================================
	# Purpose: Check for the existence of an external program.
	#
	# Details:
	#	This function is passed the name of a program or script that is s'posed 
	#	to be on the user's path, along with an option that generates a known
	#	response (the test_string).  If the passed test_string is not found in
	#	the returned output, we send a warning message and then return FALSE.
	#	The user, given a FALSE, can then cobble together a fitting response to
	#	the user.
	#
	# Example:
	#	program <- "xfm2param_nonexisting"
	#	progOptions <- "-version"
	#	test_string <- "mni_autoreg"
	#	result <- rminc.checkForExternalProgram(program, test_string, progOptions)
	#	if ( !result ) { ... }
	#
	# Note: Nothing of note, really.
	#
	# =============================================================================
	#

	# create string to submit to shell and then run it
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



rminc.readLinearXfmFile <- function(xfmFilename) {
	# =============================================================================
	# Purpose: Read the contents of a linear XFM file.
	#
	# Details:
	#	This is done by spawning the xfm2param program from the  
	#	mni_autoreg package.  As such, we need to make sure that this
	#	MNI package is installed and xfm2param is on the user's PATH.
	#
	#
	# Note: Nothing of note, really.
	#
	# =============================================================================
	#

	# check for the external program
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








