#
# =============================================================================
# Purpose: MincSlice Class Definition
#
# Description:
#	OK, here is the conceptual overview.  Volumes are read in their entirety 
#	using read volume, placing the entire 3D volume into an array.  Slices 
#	within this array are then read from and written to this 3D array
#	using the methods defined in this class. Sound reasonable?
#
#	Also, note the semantics: "read/write" prefix indicates *file* IO,
#	whereas a "get/put" prefix suggests manipulation within *memory*.
#
# TODO: 
#
#
# Notes: 
#	A.	I experimented with the idea of inheriting from class "matrix".  Doing
#		this allows you to initialize the data part of the object simply by
#		putting the data vector/matrix/array as the 2nd element in the 
#		constructor.  Doing this creates a ".Data" slot, which permits one to
#		directly access the object's data by via the usual index "[]" operators.
#		Other important methods: (1) getDataPart() gets the ".Data" data, 
#		(2) setDataPart() replaces the ".Data" slot contents.
#
# =============================================================================
#
setClass("MincSlice", 
		representation( mincInfo="MincInfo",
						sliceNumber="numeric",
						orientation="character",
						sliceIntensityRange="numeric",
						volumeType="character",
						colorMap="character",
						aspectRatio="numeric" ),
		
		prototype( sliceNumber=0,
					orientation="undefined",
					sliceIntensityRange=c(0,0),
		 			volumeType="undefined",
					colorMap="undefined",
					aspectRatio=1 ),
		contains="matrix"
)



# =============================================================================
# Methods to get/set a desired property value from a MincSlice object
#
# =============================================================================

# METHOD: mincIO.getProperty(MincSlice)
# PURPOSE: get a property from a MincSlice object
setMethod(
	"mincIO.getProperty", 
	signature=signature(mincIOobj="MincSlice", propertyId="character"),
	definition=function(mincIOobj, propertyId) {

		if ( R_DEBUG_mincIO ) cat(">> MincSlice::mincIO.getProperty() ... \n")

		# does the property exist in the MincVolumeIO object?
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

		if ( R_DEBUG_mincIO ) cat("<< MincSlice::mincIO.getProperty() ... \n")

		# return
		return(value)
	}
)



# METHOD: mincIO.setProperty(MincSlice)
# PURPOSE: set a property from a MincSlice object
setMethod(
	"mincIO.setProperty", 
	signature=signature(mincIOobj="MincSlice", propertyId="character"),
	definition=function(mincIOobj, propertyId, value) {

		if ( R_DEBUG_mincIO ) cat(">> MincSlice::mincIO.setProperty() ... \n")


		# get the variable name for the passed MincSlice object
		# ... we'll need it later to write out the updated object in place
		objName <- deparse(substitute(mincIOobj))

		# does the property exist in the MincSlice object? i.e., valid slot name passed?
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
		                      "sliceNumber",
		                      "sliceIntensityRange",
		                      "volumeType",
		                      "colorMap")


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

		if ( propertyId == "sliceNumber") mincIOobj@sliceNumber <- as.numeric(value)
		if ( propertyId == "sliceIntensityRange") mincIOobj@sliceIntensityRange <- as.numeric(value)
		if ( propertyId == "volumeType") mincIOobj@volumeType <- as.character(value)
		if ( propertyId == "colorMap") mincIOobj@colorMap <- as.character(value)


		# assign newly updated object to parent frame and then return nothing
		if ( R_DEBUG_mincIO ) {
			cat(sprintf("MincSlice::mincIO.setProperty(). Old property: %s\n", as.character(prevValue)))
			cat(sprintf("MincSlice::mincIO.setProperty(). New property: %s\n", as.character(value)))
			cat(sprintf("MincSlice::mincIO.setProperty(). Updating object: %s\n", as.character(objName)))
		}
		#
		assign(objName, mincIOobj, envir=parent.frame())
		if ( R_DEBUG_mincIO ) cat("<< MincSlice::mincIO.setProperty() ... \n")
		return(invisible())
	}
)



# =====================================================================
# Print methods
# =====================================================================

# METHOD: print(MincSlice)
# PURPOSE: print mincSlice info
setMethod(
	"print", 
	signature=signature(x="MincSlice"),
	definition=function(x) {
		if ( R_DEBUG_mincIO ) cat("MincSlice print() method ...\n")
		# assume a MincInfo object has been passed
		mincIO.printMincInfo(x@mincInfo)
	
		# display a little something about the volume data itself
		cat("\n---- Slice Specific Information ----\n")
		cat(paste("Slice number: ", x@sliceNumber, "\n", sep=""))
		cat(paste("Slice orientation: ", x@orientation, "\n", sep=""))
		xmin = x@sliceIntensityRange[1]
		xmax = x@sliceIntensityRange[2]
		cat("Min/max data values in slice:", xmin, "/", xmax, "\n")

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", x@volumeType))
		cat(sprintf("Colormap used for display: %s\n", x@colorMap))
		cat(sprintf("Aspect ratio: %s\n\n", x@aspectRatio))
	}
)


# METHOD: show(MincSlice)
# PURPOSE: print mincSlice info by simply typing "obj" at the prompt
setMethod(
	"show", 
	signature=signature(object="MincSlice"),
	definition=function(object) {
		if ( R_DEBUG_mincIO ) cat("MincSlice show() method ...\n")
		# assume a MincInfo object has been passed
		mincIO.printMincInfo(object@mincInfo)
	
		# display a little something about the volume data itself
		cat("\n---- Slice Specific Information ----\n")
		cat(paste("Slice number: ", object@sliceNumber, "\n", sep=""))
		cat(paste("Slice orientation: ", object@orientation, "\n", sep=""))
		xmin = object@sliceIntensityRange[1]
		xmax = object@sliceIntensityRange[2]
		cat("Slice real data range:", 
			sprintf("%9.3f", object@sliceIntensityRange[1]), "/", 
			sprintf("%9.3f", object@sliceIntensityRange[2]), "\n")

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", object@volumeType))
		cat(sprintf("Colormap used for display: %s\n", object@colorMap))
		cat(sprintf("Aspect ratio: %s\n\n", object@aspectRatio))
	}
)


# METHOD: plot(MincSlice)
# PURPOSE: display a minc slice
setMethod(
	"plot", 
	signature=signature(x="MincSlice", y="ANY"),
	definition=function(x, y, ...) {

		if ( R_DEBUG_mincIO ) cat("MincSlice plot() method ...\n")
		Zdim <- 1
		Ydim <- 2
		Xdim <- 3
		# do the slice plot using lattice.  
		#
		if ( x@orientation == "xSlice") {
			myPlot <- mincIO.plotSlicePretty(getDataPart(x), 
										xAxisLabel="Y Axis", 
										yAxisLabel="Z Axis", 
										aspectRatio=x@aspectRatio,
										colorMap=x@colorMap)
			print(myPlot)
		}

		if ( x@orientation == "ySlice") {
			myPlot <- mincIO.plotSlicePretty(getDataPart(x), 
										xAxisLabel="X Axis", 
										yAxisLabel="Z Axis", 
										aspectRatio=x@aspectRatio,
										colorMap=x@colorMap)
			print(myPlot)
		}

		if ( x@orientation == "zSlice") {
			myPlot <- mincIO.plotSlicePretty(getDataPart(x), 
										xAxisLabel="X Axis", 
										yAxisLabel="Y Axis", 
										aspectRatio=x@aspectRatio,
										colorMap=x@colorMap)
			print(myPlot)
		}
	}
)



mincIO.plotSlicePretty <- function(slice, xAxisLabel, yAxisLabel, aspectRatio, colorMap) {
	# =============================================================================
	# Purpose:	Common Slice Plotting Function
	#
	# Notes:	
	#	(A)		We're using Lattice, because it has fewer calories, yet tastes great.
	#	(B)		I would have preferred to draw the grid lines with panel.grid, but 
	#			sadly it's stoopid (TM). Spcifically, it finds it incomprehensible 
	#			that anyone might want to draw the grid lines where the tickmarks are.
	#
	# =============================================================================
	# 

	if ( R_DEBUG_mincIO ) cat(">> mincIO.plotSlicePretty ... \n")

	# init the colormap to use for display
	# ... first ensure that a valid colormap was specified
	clrmaps <- c("gray", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors")
	if ( !(colorMap %in% clrmaps) ) {
		warning(sprintf("Invalid colormap specified [%s].  Using grayscale colormap instead.\n", colorMap))
		colorMap <- "gray"
	}

	# generate a 256 level color gradient
	if ( colorMap == "gray") {
		colorMap.v <- eval(call(colorMap,quote(0:255/255)))
	} else {
		colorMap.v <- eval(call(colorMap,quote(256)))
	}

	# create the basic plot
	myPlot <- levelplot(slice, 
						col.regions=colorMap.v, 
						cuts=255, 
						aspect=aspectRatio,
						xlab=xAxisLabel,
						ylab=yAxisLabel,
						scales=list( x=list(tick.number=10), y=list(tick.number=10) )
						)

	# OK, we have the basic plot, now let's spice it up a bit
	myPlot <- update(myPlot,
				panel=function(...) {
					# re-do the level plot
					panel.levelplot(...)
					# add the grid to match the axis tickmarks
					# ... yes, I know, panel.grid will *not* do this correctly
					ref.line <- trellis.par.get("reference.line")
	    				panel.abline(
						v = pretty(current.panel.limits()$xlim, n = 10),
						h = pretty(current.panel.limits()$ylim, n = 10),
						col = ref.line$col,
						lty = ref.line$lty,
						lwd = ref.line$lwd
					)
				}
			)
	# done. Send it back.
	if ( R_DEBUG_mincIO ) cat("<< mincIO.plotSlicePretty ... \n")
	return(myPlot)
}



# =====================================================================
# Methods to get and return a slice from the minc volume array,
#    given slice number and orientation
#
# Input: 	MincVolumeIO object
#			slice number:	integer
#
# Output: MincSlice object
#
# Note: 
#	A.	Although it would look cleaner to have each "getSlice" method
#		call a common routine to do the common slice processing, this would
#		potentially result in the slice object needing to be copied again.
#		Thus, although less pretty, this is likely more efficient.
# =====================================================================


# METHOD: mincIO.getSliceX(MincVolumeIO, numeric)
# PURPOSE: get and return a sagittal slice from the minc volume array
setMethod(
	"mincIO.getSliceX", 
	signature=signature(mincVolume="MincVolumeIO"),
	definition=function(mincVolume, sliceNo) {
		#
		# OK, get the slice
		if ( R_DEBUG_mincIO ) cat("MincSlice method: getting X-slice from array ... \n")

		# set some useful variables
		Zdim <- 1
		Ydim <- 2
		Xdim <- 3

		# as this is a sagittal slice, compute the aspect ratio as Z-with/Y-width
		aspectRatio <- (( mincVolume@mincInfo@dimInfo$sizes[Zdim] * mincVolume@mincInfo@dimInfo$steps[Zdim] )  /
						( mincVolume@mincInfo@dimInfo$sizes[Ydim] * mincVolume@mincInfo@dimInfo$steps[Ydim] ))
		aspectRatio <- abs(aspectRatio)
#		cat(paste("... sagittal aspect ratio is ", aspectRatio, "\n"))

		# create the MincSlice object and set fields
		mincSlice <- new("MincSlice", 
							mincVolume[sliceNo,,],					# load slice data into .Data
							mincInfo=mincVolume@mincInfo, 
							sliceNumber=sliceNo,
							orientation="xSlice",
							sliceIntensityRange=range(mincVolume[sliceNo,,]),
							volumeType=mincVolume@volumeType,
							colorMap=mincVolume@colorMap,
							aspectRatio=aspectRatio)
		
		# DONE. Return the new volume array object.
		return(mincSlice)
	}
)


# METHOD: mincIO.getSliceY(MincVolumeIO, numeric)
# PURPOSE: get and return a coronal slice from the minc volume array
setMethod(
	"mincIO.getSliceY", 
	signature=signature(mincVolume="MincVolumeIO"),
	definition=function(mincVolume, sliceNo) {
		#
		# OK, get the slice
		if ( R_DEBUG_mincIO ) cat("MincSlice method: getting Y-slice from array ... \n")

		# set some useful variables
		Zdim <- 1
		Ydim <- 2
		Xdim <- 3

		# as this is a coronal slice, compute the aspect ratio as Z-with/X-width
		aspectRatio <- (( mincVolume@mincInfo@dimInfo$sizes[Zdim] * mincVolume@mincInfo@dimInfo$steps[Zdim] )  /
						( mincVolume@mincInfo@dimInfo$sizes[Xdim] * mincVolume@mincInfo@dimInfo$steps[Xdim] ))
		aspectRatio <- abs(aspectRatio)
#		cat(paste("... coronal aspect ratio is ", aspectRatio, "\n"))

		# create the MincSlice object and set fields
		mincSlice <- new("MincSlice", 
							mincVolume[,sliceNo,], 
							mincInfo=mincVolume@mincInfo, 
							sliceNumber=sliceNo,
							orientation="ySlice",
							sliceIntensityRange=range(mincVolume[,sliceNo,]),
							volumeType=mincVolume@volumeType,
							colorMap=mincVolume@colorMap,
							aspectRatio=aspectRatio)
		
		# DONE. Return the new volume array object.
		return(mincSlice)
	}
)


# METHOD: mincIO.getSliceZ(MincVolumeIO, numeric)
# PURPOSE: get and return an axial slice from the minc volume array
setMethod(
	"mincIO.getSliceZ", 
	signature=signature(mincVolume="MincVolumeIO"),
	definition=function(mincVolume, sliceNo) {
		#
		# OK, get the slice
		if ( R_DEBUG_mincIO ) cat("MincSlice method: getting Z-slice from array ... \n")

		# set some useful variables
		Zdim <- 1
		Ydim <- 2
		Xdim <- 3

		# as this is an axial slice, compute the aspect ratio as Y-with/X-width
		aspectRatio <- (( mincVolume@mincInfo@dimInfo$sizes[Ydim] * mincVolume@mincInfo@dimInfo$steps[Ydim] )  /
						( mincVolume@mincInfo@dimInfo$sizes[Xdim] * mincVolume@mincInfo@dimInfo$steps[Xdim] ))
		aspectRatio <- abs(aspectRatio)
#		cat(paste("... axial aspect ratio is ", aspectRatio, "\n"))

		# create the MincSlice object and set fields
		mincSlice <- new("MincSlice", 
							mincVolume[,,sliceNo], 
							mincInfo=mincVolume@mincInfo, 
							sliceNumber=sliceNo,
							orientation="zSlice",
							sliceIntensityRange=range(mincVolume[,,sliceNo]),
							volumeType=mincVolume@volumeType,
							colorMap=mincVolume@colorMap,
							aspectRatio=aspectRatio)
		
		# DONE. Return the new volume array object.
		return(mincSlice)
	}
)



# =====================================================================
# Methods to move an updated slice to the minc volume array
#
# Input: 	MincSlice		object
#			MincVolumeIO	object
#			slice number:	integer
#
# Output: Nutt'in
#
# Note: 
#	A.	As each slice object "knows" it's original slice number, calling this
#		method without the slice number will result in it being written
#		to its original location within the volume.  I suspect that this will
#		be the most common situation -- although I am still permitting the slice
#		to be written to a different slice within the volume.
# =====================================================================


# METHOD: mincIO.putSlice(MincSlice, MincVolumeIO, numeric)
# PURPOSE: move a slice into a volume
setMethod(
	"mincIO.putSlice", 
	signature=signature(mincSlice="MincSlice", mincVolume="MincVolumeIO"),
	definition=function(mincSlice, mincVolume, sliceNo) {
		#
		# Debug ...
		if ( R_DEBUG_mincIO ) cat("MincSlice method: writing slice back into volume array ... \n")
		Zdim <- 1
		Ydim <- 2
		Xdim <- 3


		# set the slice number to the internal value or the passed value (if valid)
		sliceNumber <- mincSlice@sliceNumber

		# if the slice number has been specified by the user, check the range
		if ( hasArg(sliceNo ) ) {
			sliceNumber <- sliceNo
			# xSpace slice?
			if (  mincSlice@orientation == "xSlice" && !(sliceNumber %in% 1:mincVolume@mincInfo@dimInfo$sizes[Xdim]) ) {
				errmsg <- sprintf("Error: Attempt to write %s slice outside of the range of the volume\n", mincSlice@orientation)
				cat(errmsg)
				errmsg <- sprintf("..... slice number: %d,  volume range for %s is [1:%d]\n", 
											sliceNumber, 
											mincSlice@orientation,
											mincVolume@mincInfo@dimInfo$sizes[Xdim])
				stop(errmsg)
			}
			# ySpace slice?
			if (  mincSlice@orientation == "ySlice" && !(sliceNumber %in% 1:mincVolume@mincInfo@dimInfo$sizes[Ydim]) ) {
				errmsg <- sprintf("Error: Attempt to write %s slice outside of the range of the volume\n", mincSlice@orientation)
				cat(errmsg)
				errmsg <- sprintf("..... slice number: %d,  volume range for %s is [1:%d]\n", 
											sliceNumber, 
											mincSlice@orientation,
											mincVolume@mincInfo@dimInfo$sizes[Ydim])
				stop(errmsg)
			}
			# zSpace slice?
			if (  mincSlice@orientation == "zSlice" && !(sliceNumber %in% 1:mincVolume@mincInfo@dimInfo$sizes[Zdim]) ) {
				errmsg <- sprintf("Error: Attempt to write %s slice outside of the range of the volume\n", mincSlice@orientation)
				cat(errmsg)
				errmsg <- sprintf("..... slice number: %d,  volume range for %s is [1:%d]\n", 
											sliceNumber, 
											mincSlice@orientation,
											mincVolume@mincInfo@dimInfo$sizes[Zdim])
				stop(errmsg)
			}
		}

		# OK, we have a valid slice number, let's write the slice back
		#
		# xSpace?
		if (  mincSlice@orientation == "xSlice" ) {
			if ( R_DEBUG_mincIO ) cat(sprintf("Writing slice %s [%d] back to volume ...\n", mincSlice@orientation, sliceNumber))
			mincVolume[sliceNumber,,] <- getDataPart(mincSlice)
		}
		#
		# ySpace?
		if (  mincSlice@orientation == "ySlice" ) {
			if ( R_DEBUG_mincIO ) cat(sprintf("Writing slice %s [%d] back to volume ...\n", mincSlice@orientation, sliceNumber))
			mincVolume[,sliceNumber,] <- getDataPart(mincSlice)
		}
		#
		# zSpace?
		if (  mincSlice@orientation == "zSlice" ) {
			if ( R_DEBUG_mincIO ) cat(sprintf("Writing slice %s [%d] back to volume ...\n", mincSlice@orientation, sliceNumber))
			mincVolume[,,sliceNumber] <- getDataPart(mincSlice)
		}

		
		# DONE. Return.
		return(mincVolume)
	}
)



# =====================================================================
# Methods to create a new empty axial slice, given a MincVolumeIO object
#
# Input: 	MincVolumeIO object
#
# Output: MincSlice object
#
# Note: None, really.
# =====================================================================


# METHOD: mincIO.makeNewSliceZ(MincVolumeIO, numeric)
# PURPOSE: create a new MincSlice object, initialized with the
#          contents of the initialization vector or zeros
setMethod(
	"mincIO.makeNewSliceZ", 
	signature=signature(mincVolume="MincVolumeIO"),
	definition=function(mincVolume, initVector) {
		#
		# OK, get the slice
		if ( R_DEBUG_mincIO ) cat("MincSlice method: creating a new Z-slice from array ... \n")

		# create a new empty z-slice (use slice 1 as a template)
		zSlice <- mincIO.getSliceZ(mincVolume, 1)

		# set the slice to zero
		sliceSize <- dim(zSlice)[1] * dim(zSlice)[2]
		
		# 
		if ( missing(initVector) ) {
			# initialization vector is not specified, so init to zeros
			zSlice <- setDataPart(zSlice, matrix(rep(0, sliceSize), nrow=dim(zSlice)[1]))
			zSlice@sliceIntensityRange <- c(0,0)
			
		} else {
			# we got an initialization vector, so validate and then use it
			# validate: vector?
			if ( class(initVector) != "numeric" ) {
				stop("Slice initialization vector must be, ummmm, ... a vector.\n")
			}
			# validate: vector of correct length?
			if ( length(initVector) != sliceSize ) {
				cat(sprintf("Slice initialization vector must contain the same number of elements as the slice\n"))
				stop(sprintf("Initialization vector: %d    Slice: %d\n", length(initVector), sliceSize))
			}
			# good. validated.  move the data in.
			zSlice <- setDataPart(zSlice, matrix(initVector, nrow=dim(zSlice)[1]))
			zSlice@sliceIntensityRange <- range(initVector)
		}
		
		# DONE. Return the new volume array object.
		return(zSlice)
	}
)








