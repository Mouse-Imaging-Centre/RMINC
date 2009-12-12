
# test to see whether files exist and are readable
isReadable <- function(filenames) {
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
	rValue <- FALSE
	if ( !file.exists(filename) ) { return(rValue) }
	if ( isMinc1(filename) ) rValue <- TRUE
	if ( isMinc2(filename) ) rValue <- TRUE
	#
	return(rValue)
}


# test to see whether a given file is minc1
isMinc1 <- function(filename) {
	rValue <- FALSE
	sysCmd <- paste("file", filename)
#	print(sysCmd)
	rtnString <- system(sysCmd, intern=TRUE)
#	print(rtnString)
	if ( grepl("NetCDF", rtnString, fixed=TRUE) ) rValue <- TRUE
	return(rValue)
}


# test to see whether a given file is minc2
isMinc2 <- function(filename) {
	rValue <- FALSE
	sysCmd <- paste("file", filename)
#	print(sysCmd)
	rtnString <- system(sysCmd, intern=TRUE)
#	print(rtnString)
	if ( grepl("Hierarchical Data Format", rtnString, fixed=TRUE) ) rValue <- TRUE
	return(rValue)
}


# convert minc1 volume to minc2
asMinc2 <- function(filename, keepName=TRUE) {
	
	# is it already minc2? Just return the input filename.
	if ( isMinc2(filename) ) return(filename)
	
	# if it isn't minc1, tell 'em and run away
	if ( !isMinc(filename) ) {
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


convertVoxelToWorld <- function(filename, voxCoords) {
	#
	
	if ( R_DEBUG_mincIO ) cat(sprintf(">>convertVoxelToWorld\n"))

	# the C routines want the coordinates 0-relative, AND in volume order
	# ... so convert first
	voxCoords <- rev(voxCoords) -1
	
	# dunno why, but the voxel coordinates are passed as doubles
	output <- .Call("convert_voxel_to_world_mincIO",
               as.character(filename),
               as.double(voxCoords), PACKAGE="RMINC")

	# return a vector of 3 doubles
	if ( R_DEBUG_mincIO ) cat(sprintf("<<convertVoxelToWorld\n"))
	return(output)
}



convertWorldToVoxel <- function(filename, worldCoords) {
	#
	# dunno why, but the voxel coordinates are passed as doubles
	if ( R_DEBUG_mincIO ) cat(sprintf(">>convertWorldToVoxel\n"))

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
	if ( R_DEBUG_mincIO ) cat(sprintf("<<convertWorldToVoxel\n"))
	return(output)
}



getDataTypes <- function() {
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



getDataClasses <- function() {
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







