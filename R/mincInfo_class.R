#
# ===================================================
# Purpose: mincInfo Class Definition
#
# Input: 
#
# Note: 
#
# ===================================================
#
setClass("MincInfo", 
		representation( volumeDataClass="numeric",
						volumeDataType="numeric",
						spaceType="character",
						nDimensions="numeric",
						dimInfo="data.frame",
						nFrames="numeric",
						timeWidths="numeric",
						timeOffsets="numeric",
						volumeIntensityRange="numeric",
						filename="character" ),
		prototype( volumeDataClass=0,
						volumeDataType=0,
						spaceType=c("init1_space", "init2_space", "init3_space"),
						nDimensions=3,
						nFrames=0,
						timeWidths=0,
						timeOffsets=0,
						volumeIntensityRange=c(0,0),
						filename="init_filename" )
)



# =====================================================================
# Print methods
# =====================================================================

# METHOD: print(MincInfo)
# PURPOSE: print MincInfo object info
setMethod("print", "MincInfo",
			function(x) {
				# assume a MincInfo object has been passed
				mincIO.printMincInfo(x)
				}
)


# METHOD: show(MincInfo)
# PURPOSE: print MincInfo object info by simply typing "obj" at the prompt
setMethod("show", "MincInfo",
			function(object) {
				# assume a MincInfo object has been passed
				mincIO.printMincInfo(object)
				}
)


# get the dimension sizes of a particular file.
mincIO.readMincInfo <- function(filename) {
	# ===================================================
	# Purpose: Read a bunch of stuff from the volume header
	#          and return a MincInfo object.  
	#
	# Input: minc2 filename
	#
	# Output: MincInfo object
	#
	# Note: 
	#
	# ===================================================
	#

	if ( R_DEBUG_mincIO ) cat(">> MincInfo::mincIO.readMincInfo() ... \n")

	# make sure that any shell wild characters are expanded
	filename <- path.expand(filename)

	# let's make sure that we have a minc2 volume (else convert)
	filename <- rminc.asMinc2(filename)

	#cat("MincInfo::readMincInfo(). About to enter the C code\n")
	if ( R_DEBUG_mincIO ) cat("MincInfo::mincIO.readMincInfo(). Entering C code\n")
	volInfo <- .Call("get_volume_info",
              as.character(filename),PACKAGE="RMINC")
	if ( R_DEBUG_mincIO ) cat("MincInfo::mincIO.readMincInfo(). Returning from C code\n")
	# print(volInfo)
	# str(volInfo)

	# create the MincInfo object and init the fields with the
	# info returned above
	mincInfo <- new("MincInfo")
	mincInfo@volumeDataClass <- volInfo$volumeDataClass
	mincInfo@volumeDataType <- volInfo$volumeDataType
	mincInfo@spaceType <- volInfo$spaceType
	mincInfo@nDimensions <- volInfo$nDimensions
	mincInfo@nFrames <- volInfo$nFrames
	mincInfo@dimInfo <- data.frame(sizes=volInfo$dimSizes, 
									steps=volInfo$dimSteps, 
									starts=volInfo$dimStarts, 
									row.names=volInfo$dimNames, 
									units=volInfo$dimUnits, stringsAsFactors=FALSE)


	# add some derived information
	# ... if we have a 4-D volume, then nFrames should have been initialized
	
	# mincInfo@nFrames <- 0
	# if ( mincInfo@nDimensions > 3 )  mincInfo@nFrames <- volInfo$nFrames


	# add the time-related info (if necessary)
	if ( mincInfo@nDimensions > 3) {
		#
		# Now let's deal with the time/timeWidths in a cludgey fashion.
		# That is, since the minc2 API is currently not able to return this info, 
		# we are going to have to use an external mincinfo call to return it.
		# Specifically, for example, although the the "time" dimension is 
		# irregularly sampled, minc2 API assumes regularly samepled and
		# happily provides the incorrect values. Neat, eh?
		#
		# The minc2 API really ought to be fixed ... but who as the time?
		# This fix, while a stinkin, rotten cludge, still works.
		# Who you lookin at? It gets the job done. Yeah, that's right.
		#
		# insert the time offsets
		cmd <- paste("mincinfo -varvalues time ", filename)
		mincInfo@timeOffsets <- as.numeric(system(cmd, intern=TRUE))
		#
		# insert the time-width information
		cmd <- paste("mincinfo -varvalues time-width ", filename)
		mincInfo@timeWidths <- as.numeric(system(cmd, intern=TRUE))
		
	}

	# store a fully-qualified path & filename
	base_filename <- basename(filename)
	dir_filename <- dirname(filename)
	#
	# have we been given a filename in the cwd without qualifier?
	if ( dir_filename == "." ) {
		filename <- file.path(getwd(), filename)
#		cat(paste("trace 1: filename = ", filename, "\n"))
	}
	# save the name away
	mincInfo@filename <- filename

	#
	if ( R_DEBUG_mincIO ) cat("<< MincInfo::mincIO.readMincInfo() ... \n")
	return(mincInfo)
}




mincIO.printMincInfo <- function(mincInfo) {
	# ===================================================
	# Purpose: Create a very user-friendly display of a
	#          MincInfo object.  
	#
	# Input: MincInfo object
	#
	# Note: Called by a number of show() methods
	#
	# ===================================================
	#

	# get the valid minc2 data types and classes
	dataClass.df <- rminc.getDataClasses()
	dataType.df <- rminc.getDataTypes()

	# great. Now print it real pretty like
	cat("\n---- Volume Specific Information ----\n")
	cat(paste("File:", mincInfo@filename, "\n"))

	# get the data class string, given the enum numeric code (and print it)
	enumCode <- which(dataClass.df$numCode == mincInfo@volumeDataClass)
	cat(paste("Interpreted data class:", dataClass.df$string[enumCode], "\n"))

	# get the data type string, given the enum numeric code (and print it)
	enumCode <- which(dataType.df$numCode == mincInfo@volumeDataType)
	cat(paste("Internal data type:", dataType.df$string[enumCode], "\n"))

	# display volume real data range
	cat(paste("Volume real data range:", 
								sprintf("%9.3f", mincInfo@volumeIntensityRange[1]), "/", 
								sprintf("%9.3f", mincInfo@volumeIntensityRange[2]), "\n"))
	

	# print the dimension-related information
	dimNames <- paste(mincInfo@dimInfo$names, collapse="  ")
	cat(paste("\nImage dimensions:", dimNames, "\n"))
	print(mincInfo@dimInfo, row.names=TRUE)

	# print the time-related information (if we have some)
	if ( mincInfo@nFrames > 0 ) {
#		cat(paste("time offsets:", mincInfo@timeOffsets, "\n"))
#		cat(paste("time widths:", mincInfo@timeWidths, "\n"))
		cat("\n----- Dynamic (Time) Information -----\n")
		tmpMat <- matrix(c(mincInfo@timeOffsets, mincInfo@timeWidths), nrow=2, byrow=TRUE, dimnames = list(c("Offsets", "Widths"), c()) )
		print(tmpMat, row.names=TRUE)
		cat("\n")
	}

#	return
}








