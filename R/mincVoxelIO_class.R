#
# =============================================================================
# Purpose: MincVoxelIO Class Definition
#
# Description:
#	Here we define that class which does IO at the smallest possible level
#	of granularity: the voxel.  As is the case with the SliceIO class, the
#	the VoxelIO class tries to make it possible to process at great number
#	of volumes (3D and 4D), without blowing out our memory.
#
#	So, it a nutshell, this class reads the real value for a specific 
#	voxel. Furthermore, the user can specify multiple 4D volumes, thus
#	allowing us to obtain the value of a specific voxel across all frames,
#	over all volumes.
#
#	This function was inspired by a similar function by Jason Lerch, buried 
#   in the very bowels of the RMINC package.
#
# ToDo: Nuttin.
#
# Notes: 
#
# =============================================================================
#
setClass("MincVoxelIO", 
		representation( mincInfo="MincInfo",
						voxelCoords="numeric",
						worldCoords="numeric",
						filenames="character" ),
		prototype( voxelCoords=c(0,0,0), 
					worldCoords=c(0,0,0),
					filenames="undefined" ),
		contains="matrix"
)




# =============================================================================
# Methods to get/set a desired property value from a MincVoxelIO object
#
# =============================================================================

# METHOD: mincIO.getProperty(MincVoxelIO)
# PURPOSE: get a property from a MincVoxelIO object
setMethod(
	"mincIO.getProperty", 
	signature=signature(mincIOobj="MincVoxelIO", propertyId="character"),
	definition=function(mincIOobj, propertyId) {

		if ( R_DEBUG_mincIO ) cat(">> MincVoxelIO::mincIO.getProperty() ... \n")

		# does the property exist in the MincVoxelIO object?
		propertyHit <- FALSE
		if ( propertyId %in% slotNames(mincIOobj) ) {
			value <- slot(mincIOobj, propertyId)
			propertyHit <- TRUE
		}

		# not found yet; does the property exist in the embedded MincInfo object?
		if ( !propertyHit ) {
			value <- mincIO.getProperty(mincIOobj@mincInfo, propertyId)
			if ( !is.null(value) ) {
				propertyHit <- TRUE
			}  
		}

		# property not found: warn and return NULL
		if ( !propertyHit ) {
			warning(sprintf("Property \"%s\" was not found in object \"%s\" [class %s]\n", 
		                              propertyId, 
		                              deparse(substitute(mincIOobj)),
		                              class(mincIOobj)))
			value <- NULL
		}

		if ( R_DEBUG_mincIO ) cat("<< MincVoxelIO::mincIO.getProperty() ... \n")

		# return
		return(value)
	}
)



# METHOD: mincIO.setProperty(MincVoxelIO)
# PURPOSE: set a property from a MincVoxelIO object
setMethod(
	"mincIO.setProperty", 
	signature=signature(mincIOobj="MincVoxelIO", propertyId="character"),
	definition=function(mincIOobj, propertyId, value) {

		if ( R_DEBUG_mincIO ) cat(">> MincVoxelIO::mincIO.setProperty() ... \n")


		# get the variable name for the passed MincVoxelIO object
		# ... we'll need it later to write out the updated object in place
		objName <- deparse(substitute(mincIOobj))

		# does the property exist in the MincVoxelIO object? i.e., valid slot name passed?
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
		valid_properties <- c("filename",
		                      "volumeIntensityRange",
		                      "filenames")

		# is the passed property settable?
		if ( !propertyId %in% valid_properties ) {
			warning(sprintf("Property \"%s\" is not settable in object \"%s\" [class %s]\n", 
		                              propertyId, 
		                              objName,
		                              class(mincIOobj)))
			return(invisible())
		}

		# set the property
		if ( propertyId == "filename") mincIOobj@mincInfo@filename <- as.character(value)
		if ( propertyId == "volumeIntensityRange") mincIOobj@mincInfo@volumeIntensityRange <- as.numeric(value)
		if ( propertyId == "volumeIntensityRange") mincIOobj@volumeIntensityRange <- as.numeric(value)

		if ( propertyId == "filename") mincIOobj@filenames <- as.character(value)


		# assign newly updated object to parent frame and then return nothing
		if ( R_DEBUG_mincIO ) {
			cat(sprintf("MincVoxelIO::mincIO.setProperty(). Old property: %s\n", as.character(prevValue)))
			cat(sprintf("MincVoxelIO::mincIO.setProperty(). New property: %s\n", as.character(value)))
			cat(sprintf("MincVoxelIO::mincIO.setProperty(). Updating object: %s\n", as.character(objName)))
		}
		#
		assign(objName, mincIOobj, envir=parent.frame())
		if ( R_DEBUG_mincIO ) cat("<< MincVoxelIO::mincIO.setProperty() ... \n")
		return(invisible())
	}
)




# =============================================================================
# Methods for print()/show() generic functions
#
# Notes:	(1) the name of the argument *must* match that used in the
#			print() generic function (that is, 'x' in this case)
# =============================================================================
# 
# print object when using "print(obj)"
setMethod(
	"print", 
	signature=signature(x="MincVoxelIO"),
	definition=function(x, ...) {

		if ( R_DEBUG_mincIO ) cat(">> MincVoxelIO::print() ... \n")

		# assume a MincInfo object has been passed
		if ( R_DEBUG_mincIO ) cat("MincVoxelIO::print() >> mincIO.printMincInfo ... \n")
		mincIO.printMincInfo(x@mincInfo)

		# display a little something about the voxel data itself
		cat("\n---- Voxel Specific Information ----\n")
		cat(sprintf("Sampled coordinate (voxel space): [ %d, %d, %d ]\n",
													x@voxelCoords[1],
													x@voxelCoords[2],
													x@voxelCoords[3]
		))

		cat(sprintf("Sampled coordinate (world space): [ %f, %f, %f ]\n",
													x@worldCoords[1],
													x@worldCoords[2],
													x@worldCoords[3]
		))

		cat(sprintf("VoxelIO matrix dimensions: [%d, %d]\n", dim(x)[1], dim(x)[2]))
		if ( length(x@filenames) > 1 ) {
			cat(sprintf("Volume Names:\n"))
			for ( ndx in 1: length(x@filenames) ) {
				cat(sprintf("%d.  %s\n", ndx, x@filenames[ndx]))
			}
		}

		if ( R_DEBUG_mincIO ) cat("<< MincVoxelIO::print() ... \n")

	}
)



# =====================================================================
# Methods to read the real value of a given voxel over frames and
#			files (as specified)
#
# Input:
#	The input will consist of one or more filenames (3D or 4D), and a voxel
#	coordinate in x,y,z order.
#
# Notes:
#	(1) coordinates are passed in x,y,z order
#	(2) pass voxel (not stx) coords
#	(3) passed coordinate values are 1-relative (matching R)
#
# =====================================================================


# METHOD: mincIO.readByVoxel(character, numeric)
# PURPOSE: read a specific voxel value across many minn volumes
setMethod(
	"mincIO.readByVoxel", 
	signature=signature(filenames="character"),
	definition=function(filenames, voxelCoords) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.readByVoxel() ... \n")

		# are all of the files readable?
		if ( !rminc.isReadable(filenames) ) {
			stop("Unreadable files found.\n")
		}

		# are all files valid minc2 volumes?
		checkedFilenames <- filenames
		for ( ndx in 1:length(filenames) ) {
			if ( !rminc.isMinc(filenames[ndx]) ) {
				stop(sprintf("Input volume [%s] does not appear to be Minc\n"))
			}
			#
			# if minc1, then we need to convert
			checkedFilenames[ndx] <- filenames[ndx]
			if ( !rminc.isMinc1(filenames[ndx]) ) {
				checkedFilenames[ndx] <- rminc.asMinc2(filenames[ndx])
			}
		}

		# OK. Now we know that all of the filenames are valid Minc2. So now we can read.
		voxelMatrix <- mincIO.readByVoxelX(filenames, voxelCoords)


		# return the VoxelIO object
		if ( R_DEBUG_mincIO ) cat("<< mincIO.readByVoxel() ... \n")
		return(voxelMatrix)

	}
)




mincIO.readByVoxelX <- function(filenames, voxelCoords) {
# =============================================================================
# Purpose: Common xovel read routine.
#
# Description:
#	Read a real voxel value of frames and volumes.
#
# TODO: 
#    (1) TBD, but there must be something.
#
#
# Notes: 
#
# =============================================================================
#
	if ( R_DEBUG_mincIO ) cat(">> mincIO.readByVoxelX() ... \n")


	nFrames <- 0
	nFiles <- length(filenames)

	# get the mincInfo details for the first volume
	mincInfoV01 <- mincIO.readMincInfo(filenames[1])
	nFrames <- mincInfoV01@nFrames

	# get the mincInfo for all files ... 
	# ... get frame count and ensure that all counts are the same
	for ( ndx in 1:nFiles ) {
		mincInfo <- mincIO.readMincInfo(filenames[ndx])
		if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): processing volume %d\n", ndx))

		# check that all 4D volumes have a consistent time dimension
		if ( mincInfo@nDimensions > 3 ) {
			if ( nFrames > 1 && rownames(mincInfo@dimInfo)[1] != "time" ) 
				cat(sprintf("Warning: First dimension of volume [%s] not equal to 'time'\n", filenames[ndx]))
		
				if ( nFrames != mincInfo@nFrames )
				stop(sprintf("Differing number of frames found at volume %s\n", filenames[ndx]))
		}
		
		# if the sampling the same across all volumes?
		# ... check the dimension sizes
		if ( !identical(mincInfoV01@dimInfo$sizes, mincInfo@dimInfo$sizes) ) {
			stop(sprintf("Volume %s appears to have different dimension sizes when compared to volume %s",
			mincInfo@filename,
			mincInfoV01@filename))
		}
		# ... check the step sizes
		if ( !identical(mincInfoV01@dimInfo$steps, mincInfo@dimInfo$steps) ) {
			stop(sprintf("Volume %s appears to have different step sizes when compared to volume %s",
			mincInfo@filename,
			mincInfoV01@filename))
		}
		# ... check the start offsets
		if ( !identical(mincInfoV01@dimInfo$starts, mincInfo@dimInfo$starts) ) {
			stop(sprintf("Volume %s appears to have different start offsets when compared to volume %s",
			mincInfo@filename,
			mincInfoV01@filename))
		}
		if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): processing volume %d ... complete\n", ndx))
	}

	# validate the coordinates
	#
	# 3 coords were passed?
	if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): at coord validation step\n"))
	if ( length(voxelCoords) != 3 ) 
		stop(sprintf("Coordinate vector must have a length of 3, not %d\n", length(voxelCoords)))


	# all coordinates are within range for the first volume?
	# ... if we have a 4-D volume, skip the first (time) index
	if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): doing 4D adjustment stuff\n"))
	offset <- 0
	if ( mincInfoV01@nDimensions > 3 )  offset <- offset +1
	for ( ndx in 1:3 ) {
		if ( voxelCoords[ndx] > mincInfoV01@dimInfo$sizes[ndx+offset] ) {
			errmsg <- sprintf("Error: Passed coordinate [%d] not within range for dimension %s\n",
			                		voxelCoords[ndx],
									mincInfoV01@dimInfo$names[ndx +offset])
			stop(errmsg)
		}
	}

	
	# make coords z,y,x order and 0-relative for C-call
	if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): making coords C-relative\n"))
	adjVoxCoords <- rev(voxelCoords) -1
	
	# OK. Now that we passsed validation, let's get the values
	if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): calling read_voxel_from_files\n"))
	voxel.mat <- .Call("read_voxel_from_files",
                  as.character(filenames),
                  as.integer(adjVoxCoords),
                  as.integer(nFiles),
                  as.integer(nFrames), PACKAGE="RMINC")

	# create the MincVoxelIO object and set assorted fields
	if ( R_DEBUG_mincIO ) cat(sprintf("Debug mincIO.readByVoxelX(): instantiating object\n"))
	voxelObj <- new("MincVoxelIO",
	 					voxel.mat,
						mincInfo=mincInfo, 
						voxelCoords=as.integer(voxelCoords),
						worldCoords=rminc.convertVoxelToWorld(filenames[1], voxelCoords),
						filenames=filenames )

	# DONE. Return the new volume array object.
	if ( R_DEBUG_mincIO ) cat("<< mincIO.readByVoxelX() ... \n")
	return(voxelObj)
	

}







