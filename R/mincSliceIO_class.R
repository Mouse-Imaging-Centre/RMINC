#
# =============================================================================
# Purpose: MincSliceIO Class Definition
#
# History:
#	Let me tell you a little story.  The initial concept of this package was
#	to do everything at the volume level.  This was supposed to include the
#	processing of 4D (time) volumes.  Specifically, reading a 4D volume would
#	result in the allocation of a 4D array, and processing would continue as
#	it does in the 3D case.
#
#	So, that was the original idea.  Unfortunately, the first time that I
#	tried to load my first 4D volume (34 frames of average resolution), the
#	malloc failed while trying to allocate ~2 GB of heap storage.  At this point
#	it became clear to me that we needed a method to deal with IO over a
#	large number of frames or a large number of files, while minimizing the
#	memory footprint. I decided to do this by allowing for a drop in the
#	IO granularity from volume to slice.  Thus, the birth of the MincSliceIO
#	class. That's my story, and I'm sticking to it.
#
# Description:
#	The purpose of this class is to make it *possible* to deal with a lot
#	of slices/volumes in one IO operation, while minimizing memory usage.
#	The primary IO class continues to be the MincVolumeIO class, which is
#	conceptualized as being the most commonly used class, and as such, is the
#	most comfortable and functional.  At this point in time, the MincSliceIO
#	class will be limited to *read-only* operation.  This may change in the 
#	future, but for the time being, the itch that I need to scratch only
#	requires read access at the slice level, so that's what I'm coding.
#
#
# ToDo: Nuttin.
#
# Notes: 
#
# =============================================================================
#
setClass("MincSliceIO", 
		representation( mincInfo="MincInfo",
						sliceNumber="numeric",
						orientation="character",			# "zSlice"
						nSlices="numeric",					# number of slices (columns) in slice matrix
						sliceOrigin="character",				# slice over volumes or frames?
						filenames="character",				# filenames for multi-volume reads
						volumeType="character",
						colorMap="character" ),
		prototype( sliceNumber=0, 
					nSlices=0, 
					nFrames=0, 
					orientation="zSlice",
					sliceOrigin="undefined",
					filenames="undefined",
					volumeType="undefined",
					colorMap="undefined" ),
		contains="matrix"
)



# =============================================================================
# Purpose:	Methods for print()/show() generic functions
#
# Notes:	(1) the name of the argument *must* match that used in the
#			print() generic function (that is, 'x' in this case)
#			(2) Note that the "print" and "show" methods can be set to 
#				produce different displays (*Cool*)
# =============================================================================


# METHOD: print(MincSliceIO)
# PURPOSE: print MincSliceIO info
setMethod(
	"print", 
	signature=signature(x="MincSliceIO"),
	definition=function(x, ...) {

		if ( R_DEBUG_mincIO ) cat(">> MincSliceIO::print() ... \n")

		# assume a MincInfo object has been passed
		if ( R_DEBUG_mincIO ) cat("MincSliceIO::print() >> printMincInfo ... \n")
		mincIO.printMincInfo(x@mincInfo)

		# display a little something about the volume data itself
		cat("\n---- Slice Specific Information ----\n")
		cat(paste("Slice number: ", x@sliceNumber, "\n", sep=""))
		cat(paste("Slice orientation: ", x@orientation, "\n", sep=""))
		cat(paste("Number of slices in matrix: ", x@nSlices, "\n", sep=""))
		cat(paste("Slice origin: ", x@sliceOrigin, "\n", sep=""))
		cat(sprintf("SliceIO matrix dimensions: [%d, %d]\n", dim(x)[1], dim(x)[2]))
		if ( length(x@filenames) > 1 ) {
			cat(sprintf("Volume Names:\n"))
			for ( ndx in 1: length(x@filenames) ) {
				cat(sprintf("%d.  %s\n", ndx, x@filenames[ndx]))
			}
		}

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", x@volumeType))
		cat(sprintf("Colormap used for display: %s\n\n", x@colorMap))

		if ( R_DEBUG_mincIO ) cat("<< MincSliceIO::print() ... \n")

	}
)


# METHOD: show(MincSliceIO)
# PURPOSE: print MincSliceIO info by simply typing "obj" at the prompt
setMethod(
	"show", 
	signature=signature(object="MincSliceIO"),
	definition=function(object) {

		if ( R_DEBUG_mincIO ) cat(">> MincSliceIO::show() ... \n")

		# assume a MincInfo object has been passed
		if ( R_DEBUG_mincIO ) cat("MincSliceIO::show() >> printMincInfo ... \n")
		mincIO.printMincInfo(object@mincInfo)

		# display a little something about the volume data itself
		cat("\n---- Slice Specific Information ----\n")
		cat(paste("Slice number: ", object@sliceNumber, "\n", sep=""))
		cat(paste("Slice orientation: ", object@orientation, "\n", sep=""))
		cat(paste("Number of slices in matrix: ", object@nSlices, "\n", sep=""))
		cat(paste("Slice origin: ", object@sliceOrigin, "\n", sep=""))
		cat(sprintf("SliceIO matrix dimensions: [%d, %d]\n", dim(object)[1], dim(object)[2]))
		if ( length(object@filenames) > 1 ) {
			cat(sprintf("Volume Names:\n"))
			for ( ndx in 1: length(object@filenames) ) {
				cat(sprintf("%d.  %s\n", ndx, object@filenames[ndx]))
			}
		}

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", object@volumeType))
		cat(sprintf("Colormap used for display: %s\n\n", object@colorMap))

		if ( R_DEBUG_mincIO ) cat("<< MincSliceIO::show() ... \n")
	
	}
)



# =====================================================================
# Methods to read a given slice across many volumes or frames
#
# Input:
#	The input will consist of either one filename, or multiple filenames.  The
#	processing is contigent on this, in the following way ...
#
#	(1) Single filename:
#		In this case, we shall assume that the user wants to extract all of the
#		frames, for the given slice, for the specified volume.  An object is
#		returned that contains an NxM matrix, where N=slice voxels, and 
#		M=slice over frames.
#
#	(2) Multiple filenames:
#		Here we assume that the user wants to retrieve the specified slice across
#		a range of 3D volumes.  The returned object contains an NxM matrix in which
#		N=slice voxels, and M=slice over volumes.
#
#
# =====================================================================


# METHOD: mincIO.readBySlice(character, numeric)
# PURPOSE: read a given slice across many volumes
setMethod(
	"mincIO.readBySlice", 
	signature=signature(filenames="character", sliceNumber="numeric"),
	definition=function(filenames, sliceNumber, ..., volumeType, colorMap) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.readBySlice() ... \n")

		# figure out whether we have a 3D or 4D volume
		if ( !rminc.isMinc(filenames[1]) ) {
			stop(sprintf("Specified filename [%s] either does not exist or is not minc"))
		}
		filename <- rminc.asMinc2(filenames[1])
		mincInfo <- mincIO.readMincInfo(filename[1])

		# set the display properties to useful defaults (if unset)
		if ( !hasArg(volumeType) ) {
			volumeType <- "anatomical"
			if ( mincInfo@nFrames > 0 )  volumeType <- "functional"
		}
		if ( !hasArg(colorMap) ) {
			colorMap <- "gray"
			if ( volumeType == "functional")  colorMap <- "rainbow"
		}


		# case (1) -- read all frames for a single 4D volume
		if ( mincInfo@nDimensions == 4 ) {

			# set some useful variables (instead of hard-coding in the function body)
			tDim <- 1
			zDim <- 2
			yDim <- 3
			xDim <- 4
			
			# validate the slice number
			if ( !(sliceNumber %in% 1:mincInfo@dimInfo$sizes[zDim]) ) {
				stop(sprintf("Specified slice number [%d] does not fall with the valid range of volume %s [1:%d]",
									sliceNumber, mincInfo@filename, mincInfo@dimInfo$sizes[zDim]))
			}

			# set start indices and counts to read an entire volume, then read
			startIndices <- rep(0, mincInfo@nDimensions)
			startIndices[zDim] <- sliceNumber -1								# make the index 0-relative
			#
			dimSizes <- mincInfo@dimInfo$sizes
			dimSizes[zDim] <- 1													# only get one slice

			# do the read
			sliceMatrix <- .Call("read_hyperslab",
		               as.character(mincInfo@filename),
		               as.integer(startIndices),
		               as.integer(dimSizes),
		               as.integer(mincInfo@nDimensions), PACKAGE="RMINC")

			# make into a 2D matrix [x,y]
			dim(sliceMatrix) <- c(mincInfo@dimInfo$sizes[xDim] * mincInfo@dimInfo$sizes[yDim], mincInfo@dimInfo$sizes[tDim])

			# create the MincSliceIO object and set assorted fields
			mincSlice <- new("MincSliceIO",
			 					sliceMatrix,
								mincInfo=mincInfo, 
								sliceNumber=sliceNumber,
								filenames=mincInfo@filename,
								nSlices=mincInfo@dimInfo$sizes[tDim],
								sliceOrigin="Slices over frames",
								volumeType=volumeType,
								colorMap=colorMap)
		}


		# case (2) -- read a specific slice from a series of 3D volumes
		if ( mincInfo@nDimensions == 3 ) {

			# set some useful variables (instead of hard-coding in the function body)
			zDim <- 1
			yDim <- 2
			xDim <- 3
			
			# validate files and convert to minc2 (if necessary)
			nVolumes <- length(filenames)
			for ( ndx in 1:nVolumes ) {
				if ( !rminc.isMinc(filenames[ndx] ) ) {
					stop(sprintf("Specified filename [%s] either does not exist or is not minc"))
				}
				filenames[ndx] <- rminc.asMinc2(filenames[ndx])
			}
			
			# Initialize the MincInfo object for the first volume 
			# ... yes, we are going to assume that Mr. User is smart enough not to feed
			# ... us volumes with different sampling
			if ( R_DEBUG_mincIO ) cat(sprintf("About to read %s\n", filenames[1]))
			mincInfo <- mincIO.readMincInfo(filenames[1])
			if ( R_DEBUG_mincIO ) cat(sprintf("Reading complete\n"))
			
			# don't allow reading anything other than a 3D volume
			if ( mincInfo@nDimensions != 3) {
				cat(sprintf("Error: Sorry, but readBySlice can only read 3D volumes when specifying multiple input volumes\n"))
				stop(sprintf("Your volume [%s] has %d dimensions\n", mincInfo@filename, mincInfo@nDimensions))
			}
			# validate the slice number
			if ( !(sliceNumber %in% 1:mincInfo@dimInfo$sizes[zDim]) ) {
				stop(sprintf("Specified slice number [%d] does not fall with the valid range of volume %s [1:%d]",
										sliceNumber, mincInfo@filename, mincInfo@dimInfo$sizes[zDim]))
			}

			# create an empty MincSliceIO object and set assorted fields
			sliceLength <- mincInfo@dimInfo$sizes[xDim] * mincInfo@dimInfo$sizes[yDim]
			mincSlice <- new("MincSliceIO",
			 					matrix(rep(0, sliceLength*nVolumes), ncol=nVolumes),
								mincInfo=mincInfo, 
								sliceNumber=sliceNumber,
								filenames=filenames,
								nSlices=nVolumes,
								sliceOrigin="Slices over volumes",
								volumeType=volumeType,
								colorMap=colorMap)
			
			# set start indices and counts to read the specified zSlice
			startIndices <- rep(0, mincInfo@nDimensions)
			startIndices[zDim] <- sliceNumber -1								# make the index 0-relative
			#
			dimSizes <- mincInfo@dimInfo$sizes
			dimSizes[zDim] <- 1													# only get one slice

			# loop over all volumes and read the slice
			for ( ndx in 1:length(filenames) ) {
				if ( R_DEBUG_mincIO ) cat(sprintf("setting slice array column %d from volume %s\n", ndx, filenames[ndx]))
				mincSlice[,ndx] <- .Call("read_hyperslab",
										as.character(filenames[ndx]),
										as.integer(startIndices),
										as.integer(dimSizes),
										as.integer(mincInfo@nDimensions), PACKAGE="RMINC")
			}
		}

		# return the slice object
		if ( R_DEBUG_mincIO ) cat("<< mincIO.readBySlice() ... \n")
		return(mincSlice)

	}
)



# =====================================================================
# Methods to get and return a slice/frame from the minc slice matrix,
#         given slice number
#
# Input: 	MincSliceIO object
#			slice index:	integer index to given matrix column
#
# Output: MincSlice object
#
# Note: 
#	None.
# =====================================================================


# METHOD: mincIO.getSliceFromSliceArray(MincSliceIO, numeric)
# PURPOSE: return a specific slice from a slice array
setMethod(
	"mincIO.getSliceFromSliceArray", 
	signature=signature(mincSliceMatrix="MincSliceIO"),
	definition=function(mincSliceMatrix, sliceIndex) {
		#
		# OK, create a new slice object and init it
		if ( R_DEBUG_mincIO ) cat("MincSliceIO method: getting slice/frame from minc Slice Matrix ... \n")

		# set some useful variables (instead of hard-coding in the function body)
		offset <- mincSliceMatrix@mincInfo@nDimensions -3
		zDim <- offset +1
		yDim <- offset +2
		xDim <- offset +3

		# compute axial slice aspect ratio as X-with/Y-width
		xSize <- mincSliceMatrix@mincInfo@dimInfo$sizes[xDim] * mincSliceMatrix@mincInfo@dimInfo$steps[xDim]
		ySize <- mincSliceMatrix@mincInfo@dimInfo$sizes[yDim] * mincSliceMatrix@mincInfo@dimInfo$steps[yDim]
		aspectRatio <- ySize / xSize
		if ( R_DEBUG_mincIO ) {
			cat(sprintf("x-axis length: %d  step: %f   Total width: %f\n", 
									mincSliceMatrix@mincInfo@dimInfo$sizes[xDim],
									mincSliceMatrix@mincInfo@dimInfo$steps[xDim],
									xSize))
		}
		if ( R_DEBUG_mincIO ) {
			cat(sprintf("y-axis length: %d  step: %f   Total width: %f\n", 
									mincSliceMatrix@mincInfo@dimInfo$sizes[yDim],
									mincSliceMatrix@mincInfo@dimInfo$steps[yDim],
									ySize))
		}
		if ( R_DEBUG_mincIO ) cat(paste("... axial aspect ratio is ", aspectRatio, "\n"))
		if ( R_DEBUG_mincIO ) cat(paste("... nrow:", mincSliceMatrix@mincInfo@dimInfo$sizes[yDim], "\n"))
		if ( R_DEBUG_mincIO ) cat(paste("... ncol:", mincSliceMatrix@mincInfo@dimInfo$sizes[xDim], "\n"))

		# create the MincSlice object and set fields
		mincSlice <- new("MincSlice", 
							matrix(mincSliceMatrix[,sliceIndex], 
								nrow=mincSliceMatrix@mincInfo@dimInfo$sizes[xDim],
								ncol=mincSliceMatrix@mincInfo@dimInfo$sizes[yDim]),
							mincInfo=mincSliceMatrix@mincInfo, 
							sliceNumber=mincSliceMatrix@sliceNumber,
							orientation="zSlice",
							sliceIntensityRange=range(mincSliceMatrix[,sliceIndex]),
							volumeType=mincSliceMatrix@volumeType,
							colorMap=mincSliceMatrix@colorMap,
							aspectRatio=aspectRatio)
		
		# DONE. Return the new volume array object.
		return(mincSlice)
	}
)






