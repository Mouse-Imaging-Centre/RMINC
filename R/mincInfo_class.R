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



# =============================================================================
# Methods to get/set a desired property value from a MincInfo object
#
# =============================================================================

# METHOD: mincIO.getProperty(MincInfo)
# PURPOSE: get a property from a MincInfo object
setMethod(
	"mincIO.getProperty", 
	signature=signature(mincIOobj="MincInfo", propertyId="character"),
	definition=function(mincIOobj, propertyId) {

		if ( R_DEBUG_mincIO ) cat(">> MincInfo::mincIO.getProperty() ... \n")

		# does the property exist in the MincInfo object? i.e., valid slot name passed?
		propertyHit <- FALSE
		if ( propertyId %in% slotNames(mincIOobj) ) {
			value <- slot(mincIOobj, propertyId)
			propertyHit <- TRUE
		}

		# does the property exist in the dimInfo slot?
		if ( propertyId %in% names(slot(mincIOobj,"dimInfo")) ) {
			dimX <- slot(mincIOobj,"dimInfo")[propertyId]
			
			# info is returned as a data.frame in z,y,x order -- 
			# ... let's convert to a named vector in x,y,z order
			value <- dimX[rev(1:nrow(dimX)),]		# now it's an xyz vector
			names(value) <- rev(row.names(dimX))	# ... and now it's named
			propertyHit <- TRUE
		}  

		# property not found: warn and return NULL
		if ( !propertyHit ) {
			warning(sprintf("Property \"%s\" was not found in object \"%s\" [class %s]\n", 
		                              propertyId, 
		                              deparse(substitute(mincIOobj)),
		                              class(mincIOobj)))
			value <- NULL
		}

		if ( R_DEBUG_mincIO ) cat("<< MincInfo::mincIO.getProperty() ... \n")

		# return
		return(value)
	}
)



# METHOD: mincIO.setProperty(MincInfo)
# PURPOSE: set a property from a MincInfo object
setMethod(
	"mincIO.setProperty", 
	signature=signature(mincIOobj="MincInfo", propertyId="character"),
	definition=function(mincIOobj, propertyId, value) {

		if ( R_DEBUG_mincIO ) cat(">> MincInfo::mincIO.setProperty() ... \n")


		# get the variable name for the passed MincInfo object
		# ... we'll need it later to write out the updated object in place
		objName <- deparse(substitute(mincIOobj))

		# does the property exist in the MincInfo object? i.e., valid slot name passed?
		propertyHit <- FALSE
		prevValue <- mincIO.getProperty(mincIOobj, propertyId)
		if ( !is.null(prevValue) ) {
				propertyHit <- TRUE
		}

		# property not found: warn and return NULL
		if ( !propertyHit ) {
			warning(sprintf("Property \"%s\" was not found in object \"%s\" [class %s]\n", 
		                              propertyId, 
		                              objName,
		                              class(mincIOobj)))
			return(invisible())
		}

		# setting new property
		#
		# define settable properties (for this class)
		valid_properties <- c("filename", "volumeIntensityRange")

		# is the passed property settable?
		if ( !propertyId %in% valid_properties ) {
			warning(sprintf("Property \"%s\" is not settable in object \"%s\" [class %s]\n", 
		                              propertyId, 
		                              objName,
		                              class(mincIOobj)))
			return(invisible())
		}

		# set the property
		if ( propertyId == "filename") mincIOobj@filename <- as.character(value)
		if ( propertyId == "volumeIntensityRange") mincIOobj@volumeIntensityRange <- as.numeric(value)


		# assign newly updated object to parent frame and then return nothing
		if ( R_DEBUG_mincIO ) {
			cat(sprintf("MincInfo::mincIO.setProperty(). Old property: %s\n", as.character(prevValue)))
			cat(sprintf("MincInfo::mincIO.setProperty(). New property: %s\n", as.character(value)))
			cat(sprintf("MincInfo::mincIO.setProperty(). Updating object: %s\n", as.character(objName)))
		}
		#
		assign(objName, mincIOobj, envir=parent.frame())
		if ( R_DEBUG_mincIO ) cat("<< MincInfo::mincIO.setProperty() ... \n")
		return(invisible())
	}
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








