#
# =============================================================================
# Purpose: MincVolumeIO Class Definition
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
setClass("MincVolumeIO", 
		representation( mincInfo="MincInfo",
						volumeIntensityRange="numeric",
						frameNumber="numeric",
						volumeType="character",
						colorMap="character" ),
		prototype( volumeIntensityRange=c(0,0),
					frameNumber=0,
		 			volumeType="undefined",
					colorMap="undefined" ),
		contains="array"
)



# =============================================================================
# Purpose:	Methods for print()/show() generic functions
#
# Notes:	(1) the name of the argument *must* match that used in the
#			print() generic function (that is, 'x' in this case)
#			(2) Note that the "print" and "show" methods can be set to 
#				produce different displays (*Cool*)
# =============================================================================
# 
# print object when using "print(obj)"
setMethod(
	"print", 
	signature=signature(x="MincVolumeIO"),
	definition=function(x, ...) {

		if ( R_DEBUG_mincIO ) cat(">> MincVolumeIO::print() ... \n")

		# assume a MincInfo object has been passed
		if ( R_DEBUG_mincIO ) cat("MincVolumeIO::print() >> mincIO.printMincInfo() ... \n")
		mincIO.printMincInfo(x@mincInfo)

		# print out frame info (if we got it)
		if ( x@frameNumber > 0 ) {
			cat("\n---- Frame Specific Information ----\n")
			cat(sprintf("Frame Number: %d / %d\n", x@frameNumber, x@mincInfo@nFrames))
		}

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", x@volumeType))
		cat(sprintf("Colormap used for display: %s\n\n", x@colorMap))

		if ( R_DEBUG_mincIO ) cat("<< MincVolumeIO::print() ... \n")

	}
)

# print object when simply typing "obj" at the prompt
setMethod(
	"show", 
	signature=signature(object="MincVolumeIO"),
	definition=function(object) {

		if ( R_DEBUG_mincIO ) cat(">> MincVolumeIO::show() ... \n")

		# assume a MincInfo object has been passed
		if ( R_DEBUG_mincIO ) cat("MincVolumeIO::show() >> mincIO.printMincInfo() ... \n")
		mincIO.printMincInfo(object@mincInfo)

		# print out frame info (if we got it)
		if ( object@frameNumber > 0 ) {
			cat("\n---- Frame Specific Information ----\n")
			cat(sprintf("Frame Number: %d / %d\n", object@frameNumber, object@mincInfo@nFrames))
		}

		# the following fields mostly effect display
		cat("\n---- Volume Display-Related Properties ----\n")
		cat(sprintf("Volume type: %s\n", object@volumeType))
		cat(sprintf("Colormap used for display: %s\n\n", object@colorMap))

		if ( R_DEBUG_mincIO ) cat("<< MincVolumeIO::show() ... \n")
	
	}
)



# use the plot() generic to display the image
setMethod(
	"plot", 
	signature=signature(x="MincVolumeIO", y="ANY"),
	definition=function(x, y, ...) {

		if ( R_DEBUG_mincIO ) cat(">> MincVolumeIO::plot() ... \n")


		# set some useful variables (instead of hard-coding in the function body)
		offset <- x@mincInfo@nDimensions -3
		zDim <- offset +1
		yDim <- offset +2
		xDim <- offset +3
		nSlicesDisplay <- 21
		imgLayout <- c(7,3)

		# we are not going to plot every slice in the volume (takes too long)
		# so, lets select a number of equally-spaced slices
		#
		# don't try to display more slices than in volume
		nSlicesVolume <- x@mincInfo@dimInfo$sizes[zDim]
		nSlicesDisplay <- min(nSlicesVolume, nSlicesDisplay)
		selectedSlices <- round(seq(1, nSlicesVolume, length.out=nSlicesDisplay))
		sliceIds <- paste("slice", selectedSlices)



		# sub-sample the 3D volume, selecting only the slices for plotting
		plotVol <- x[,,selectedSlices]

		# compute the aspect ratio for axial slices
		aspectRatio <- (( x@mincInfo@dimInfo$sizes[yDim] * x@mincInfo@dimInfo$steps[yDim] )  /
						( x@mincInfo@dimInfo$sizes[xDim] * x@mincInfo@dimInfo$steps[xDim] ))
		aspectRatio <- abs(aspectRatio)


		# init the colormap to use for display
		# ... first ensure that a valid colormap was specified
		clrmaps <- c("gray", "rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors")
		if ( !(x@colorMap %in% clrmaps) ) {
			warning(sprintf("Invalid colormap specified [%s].  Using grayscale colormap instead.\n", x@colorMap))
			colorMap <- "gray"
		} else {
			colorMap <- x@colorMap
		}

		# generate a 256 level color gradient
		if ( colorMap == "gray") {
			colorMap.v <- eval(call(colorMap,quote(0:255/255)))
		} else {
			colorMap.v <- eval(call(colorMap,quote(256)))
		}
		
		# create the image plot (using lattice for the time being)
		myLevelPlot <- levelplot(plotVol, 
								layout=imgLayout,
								col.regions=colorMap.v, 
								cuts=255, 
								as.table=TRUE, 
								labels=FALSE,
								aspect=aspectRatio,
								xlab=NULL,
								ylab=NULL,
								main=list(basename(x@mincInfo@filename), col="yellow"),
								scales=list(draw=FALSE),
								between=list(x=c(0.5), y=c(0.5)),
								strip=strip.custom(style=1, 
													bg="black", 
													factor.levels = sliceIds, 
													par.strip.text=list(col="white", cex=.8)),
								colorkey=list(labels=list(col="yellow"))
						)

		# create the histogram (also with lattice)
		myHistPlot <- histogram(x[,,],
								col="yellow",
								xlab=list("Volume Real Intensity Values", col="yellow"),
								ylab=list(col="yellow"),
								nint=50,
								scales=list(x=list(tick.number=10), col="yellow", col.line="yellow")
						)

		# use grid graphics to plot the volume slices on top and the histogram below
		#
		# plot a black background
		grid.rect(gp=gpar(fill="black"))

		# plot the slices (create viewport, then plot)
		image.vp <- viewport(x=unit(.0,"npc"), 
								y=unit(0.40,"npc"), 
								width=unit(1.0,"npc"), 
								height=unit(0.60,"npc"), 
								just=c("left", "bottom"),
								name="image.vp")
		pushViewport(image.vp)
#		grid.rect()
		print(myLevelPlot, newpage=FALSE)

		# now plot the histogram (pop up one viewport level first)
		popViewport()
		histogram.vp <- viewport(x=unit(0,"npc"), 
									y=unit(0,"npc"), 
									width=unit(1.0,"npc"), 
									height=unit(0.4,"npc"), 
									just=c("left", "bottom"),
									name="histogram.vp")
		pushViewport(histogram.vp)
#		grid.rect()
		print(myHistPlot, newpage=FALSE)

		if ( R_DEBUG_mincIO ) cat("<< MincVolumeIO::plot() ... \n")
		
	}
)




# =============================================================================
# Purpose:	Create a generic function for each "readVolume" function, and then
#			add various methods that attach to that generic
#
# =============================================================================
# 
setGeneric( 
	name="mincIO.readVolume", 
	def = function(object, frameNo=0, ..., volumeType, colorMap) { standardGeneric("mincIO.readVolume") }
) 

# read the volume, by passing a MincInfo object
setMethod(
	"mincIO.readVolume", 
	signature=signature(object="MincInfo"),
	definition=function(object, frameNo, ..., volumeType, colorMap) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.readVolume() ... \n")

		#
		# we've been passed a MincInfo object, so we can start the read directly
		if ( R_DEBUG_mincIO ) cat("mincIO.readVolume() >> mincIO.readVolumeX() ... \n")
		mincVolume <- mincIO.readVolumeX(object, frameNo, volumeType, colorMap)
	}
)

# read the volume, by passing the "filename" of the volume to read
setMethod(
	"mincIO.readVolume", 
	signature=signature(object="character"),
	definition=function(object, frameNo, ..., volumeType, colorMap) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.readVolume() ... \n")

		# mincId is a filename
		filename <- object
	
		# make sure that the input file is minc/minc2
		if ( !rminc.isMinc(filename) ) {
			errmsg <- paste("File", filename, "does not appear to be a valid minc volume")
			stop(errmsg)
		}
		# convert if needed
		filename <- rminc.asMinc2(filename)

		# get some volume info
		mincInfo <- mincIO.readMincInfo(filename)
	
		# do the read
		if ( R_DEBUG_mincIO ) cat("mincIO.readVolume() >> mincIO.readVolumeX() ... \n")
		mincVolume <- mincIO.readVolumeX(mincInfo, frameNo, volumeType, colorMap)

		if ( R_DEBUG_mincIO ) cat("<< mincIO.readVolume() ... \n")
		return(mincVolume)
	}
)



mincIO.readVolumeX <- function(mincInfo, frameNo, volumeType, colorMap) {
#
# =============================================================================
# Purpose: Common 3D volume read routine.
#
# Description:
#	Read a 3D volume or a single 4D frame into a 3D array and return 
#	a MincVolumeIO object.
#
# TODO: 
#    (1) I am copying the output of read_hyperslab into the
#        MincVolumeIO object.  I should consider calling
#        a C function "read_volume" which calls read_hyperslab,
#        but returns an instantiated MincVolumeIO object.
#
#
# Notes: 
#
# =============================================================================
#
	if ( R_DEBUG_mincIO ) cat(">> mincIO.readVolumeX() ... \n")

	# If this is a 4d volume, make sure that a valid frame number has been passed
	if ( mincInfo@nDimensions == 4 ) {
		if ( !hasArg(frameNo) ) {
			stop(sprintf("Reading a 4D volume requires a valid frame number to be passed"))
		} else {
			if ( frameNo < 1 || frameNo > mincInfo@nFrames) {
				stop(sprintf("Invalid frame number [%d] passed.  Must be between [1,%d]", frameNo, mincInfo@nFrames))
			}
		}
	}


	# set the display properties to useful defaults (if unset)
	if ( !hasArg(volumeType) ) {
		volumeType <- "anatomical"
		if ( mincInfo@nFrames > 0 )  volumeType <- "functional"
	}
	if ( !hasArg(colorMap) ) {
		colorMap <- "gray"
		if ( volumeType == "functional")  colorMap <- "rainbow"
	}

	# set start indices and counts to read an entire volume, then read
	# make some adjustments in the case of 4 dimensions (time is dim #1)
	if (  mincInfo@nDimensions == 3 )  { 
		startIndices <-  rep(0, 3)
		zyxSizes <- mincInfo@dimInfo$sizes
		readSizes <- zyxSizes
	} else {
		startIndices <- c(frameNo -1, rep(0, 3))
		tzyxSizes <- mincInfo@dimInfo$sizes
		zyxSizes <- tzyxSizes[2:4]
		readSizes <- c(1, zyxSizes)
	}
	
	# do the read
	volume <- .Call("read_hyperslab",
               as.character(mincInfo@filename),
               as.integer(startIndices),
               as.integer(readSizes),
               as.integer(mincInfo@nDimensions), PACKAGE="RMINC")


	# make into an array
	# note that while minc expects the slowest varying dimension to be first,
	# R expects it to be last. So, we need to reverse the dim order here.
	dim(volume) <- rev(zyxSizes)

	# update volume-related fields in the MincInfo object before storing it
	mincInfo@volumeIntensityRange = c(min(volume), max(volume))

	# create the MincVolumeIO object and set assorted fields
	mincVolume <- new("MincVolumeIO",
	 					volume,
						mincInfo=mincInfo, 
						volumeIntensityRange=mincInfo@volumeIntensityRange,
						frameNumber=frameNo)

	# set display-related properties
	if ( hasArg(volumeType) ) {
		mincVolume@volumeType <- volumeType
	}
	if ( hasArg(colorMap) ) {
		mincVolume@colorMap <- colorMap
	}

	# DONE. Return the new volume array object.
	if ( R_DEBUG_mincIO ) cat("<< readVolumeX ... \n")
	return(mincVolume)
	
}





# =============================================================================
# Purpose:	Create a generic function for the "writeVolume" function, and then
#			add various methods that attach to that generic
#
# =============================================================================
# 
setGeneric( 
	name="mincIO.writeVolume", 
	def = function(object, filename=missing, clobber=FALSE) { standardGeneric("mincIO.writeVolume") }
) 

# write the volume, by passing a MincVolumeIO object, and the desired output filename
setMethod(
	"mincIO.writeVolume", 
	signature=signature(object="MincVolumeIO"),
	definition=function(object, filename, clobber) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.writeVolume() ... \n")

		# if no filename has been specified, take it from the object to be written
		if ( missing(filename) ) {
			filename <- object@mincInfo@filename
		}

		# if the file already exists AND clobber=TRUE, 
		# ... then delete the file, as we are going to replace it
		if ( file.exists(filename) ) {
			if ( clobber ) {
				file.remove(filename)
			} else {
				stop("Attempting to over-write pre-existing volume without Clobber flag set\n")
			}
		}
		
		# minc2 API does not like unexpanded paths (i.e., with "~", or "..", etc)
		# so let's ensure that it's expanded
		filename <- path.expand(filename)

		# OK, yes I know that the actual volume write could be inserted directly
		# in here, but at this stage, I'm not sure if I might overload this 
		# function more in the near future. So, let's just init and then call 
		# the function that does all of the real work.
		#
		# note that we really only care about the "side effect" here (i.e. writing
		# out the volume), so we do not need to capture a return value.

		if ( R_DEBUG_mincIO ) cat("mincIO.writeVolume() >> mincIO.writeVolumeX() ... \n")
			mincIO.writeVolumeX(object, filename)
		if ( R_DEBUG_mincIO ) cat("<< mincIO.writeVolume() ... \n")
		
		
	}
)




mincIO.writeVolumeX <- function(mincVolume, filename) {
#
# =============================================================================
# Purpose: Common 3D volume write routine.
#
# Description:
#	Write a 3D volume from the 3D array (i.e. the MincVolumeIO object).
#
# TODO: 
#    (1) Ummm, general efficiency issues should be explored.
#
#
#
# Notes: 
#
# =============================================================================
#
	if ( R_DEBUG_mincIO ) cat(">> mincIO.writeVolumeX() ... \n")


	# set start indices and counts to read an entire volume, then read
	if ( R_DEBUG_mincIO ) cat(">> mincIO.writeVolumeX() >> mincIO.write_volume() ... \n")
	callStatus <- .Call("write_volume",
               as.character(filename),
               as.integer(mincVolume@mincInfo@nDimensions),
               as.integer(mincVolume@mincInfo@dimInfo$sizes),
               as.double(mincVolume@mincInfo@dimInfo$starts),
               as.double(mincVolume@mincInfo@dimInfo$steps),
               as.integer(mincVolume@mincInfo@volumeDataType),
               as.double( c( min(mincVolume), max(mincVolume) ) ),
               as.double( getDataPart(mincVolume) ), PACKAGE="RMINC")

	# DONE. Return nothing.
	if ( R_DEBUG_mincIO ) cat("<< mincIO.writeVolumeX() ... \n")
	return
	
}





# =============================================================================
# Purpose:	Create a generic function for the "makeNewVolume" function, and then
#			add various methods that attach to that generic
#	Note:
#		This generic is heavily over-loaded, linking to 3 different methods.
#		Please remember the following S4 info:
#		(1) the args in the generic and the methods must match exactly
#		(2) the args in the signature are only used for method selection, so we
#			only need to specify those args that provide a unqiue signature.
#		(3) the signature for "likeTemplate" and "likeFile" are both 
#			c("character", "character"), so the args names need to be specified
#			in the signature.  I have explicitly set "missing" where appropriate
#			to emphasize the differences in these signatures.
# =============================================================================
# 
setGeneric( 
	name="mincIO.makeNewVolume", 
	def = function(filename=filename, 
					dimLengths, 
					dimSteps, 
					dimStarts,
					likeTemplate,
					likeFile) { standardGeneric("mincIO.makeNewVolume") }
) 


# Create a new, empty MincVolumeIO object, by passing a filename and the sampling details
setMethod(
	"mincIO.makeNewVolume", 
	signature=signature(filename="character", 
						dimLengths="numeric", 
						dimSteps="numeric", 
						dimStarts="numeric"),
	definition=function(filename=filename, 
						dimLengths=NULL, dimSteps=NULL, dimStarts=NULL,
						likeTemplate=NULL,
						likeFile=NULL) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.makeNewVolume() ... \n")

		# use the passed parameters to create a MincInfo object
		# ... use some reasonable defaults
		mincInfo <- new("MincInfo")

		# set filename
		mincInfo@filename <- filename

		# set the default data class to REAL
		dataClass.df <- rminc.getDataClasses()
		enumCode <- which(dataClass.df$string == "REAL")
		mincInfo@volumeDataClass <- dataClass.df$numCode[enumCode]

		# set the default data storage type to 16-bit unsigned integer
		dataType.df <- rminc.getDataTypes()
		enumCode <- which(dataType.df$code == "MI_TYPE_USHORT")
		mincInfo@volumeDataType <- dataType.df$numCode[enumCode]

		# options: "native____", "talairach_", etc
		mincInfo@spaceType <- "talairach_"

		mincInfo@nDimensions <- 3

		# set dim-related goodness
		# Note: In order to keep the user experience consistent, we are requiring 
		# dimensional information be entered in xyz order.
		# Of course, the minc conventions is to have the slowest  varying dimension to be first,
		# ... so we will have to reverse the user-specified (input) order.
		# ... That is, xyz --> zyx
		mincInfo@dimInfo <- data.frame(sizes=rev(dimLengths), 
										steps=rev(dimSteps), 
										starts=rev(dimStarts), 
										row.names=c("zspace", "yspace", "xspace"), 
										units=c("mm", "mm", "mm"), stringsAsFactors=FALSE)

		# use this mincInfo to create a new, empty volume
		mincVolume <- mincIO.makeNewVolumeX(mincInfo)

		# DONE. Return MincInfo object.
		if ( R_DEBUG_mincIO ) cat("<< mincIO.makeNewVolume() ... \n")
		return(mincVolume)

		
	}
)



# Create a new, empty MincVolumeIO object, by passing a filename and the sampling details
setMethod(
	"mincIO.makeNewVolume", 
	signature=signature(filename="character", 
						likeTemplate="character",
						likeFile="missing"),
	definition=function(filename=filename, 
						dimLengths=NULL, dimSteps=NULL, dimStarts=NULL,
						likeTemplate="icbm152",
						likeFile=NULL) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.makeNewVolume()/templates ... \n")


		# make sure that a valid template volume has been specified
		if ( likeTemplate != "icbm152" &&
		 	likeTemplate != "mni305linear" &&
		 	likeTemplate != "mni305PET" ) {
			stop(sprintf("Error: Specified template volume [%s] is not supported.\n", likeTemplate));
		}

		# Yes. Very nice.  Now create an empty MincInfo object
		mincInfo <- new("MincInfo")


		#
		# Now let's tailor the MincInfo to reflect a given known template volume
		#

		# (1) volume should look like an ICBM152-sampled volume
		if ( likeTemplate == "icbm152" ) {

			# set the MincInfo fields
			mincInfo@filename <- filename

			dataClass.df <- rminc.getDataClasses()
			enumCode <- which(dataClass.df$string == "REAL")
			mincInfo@volumeDataClass <- dataClass.df$numCode[enumCode]

			dataType.df <- rminc.getDataTypes()
			enumCode <- which(dataType.df$code == "MI_TYPE_SHORT")
			mincInfo@volumeDataType <- dataType.df$numCode[enumCode]

			mincInfo@spaceType <- "talairach_"
			mincInfo@nDimensions <- 3
			# dim info ordering is zyx
			mincInfo@dimInfo <- data.frame(sizes=c(181, 217, 181), 
											steps=c(1, 1, 1), 
											starts=c(-72, -126, -90), 
											row.names=c("zspace", "yspace", "xspace"), 
											units=c("mm", "mm", "mm"), stringsAsFactors=FALSE)
		}


		# (2) volume should look like the MNI 305 linear volume 
		#	(as found in the mni-models_average305-lin-1.0 package)
		if ( likeTemplate == "mni305linear" ) {

			# set the MincInfo fields
			mincInfo@filename <- filename

			dataClass.df <- rminc.getDataClasses()
			enumCode <- which(dataClass.df$string == "REAL")
			mincInfo@volumeDataClass <- dataClass.df$numCode[enumCode]

			dataType.df <- rminc.getDataTypes()
			enumCode <- which(dataType.df$code == "MI_TYPE_UBYTE")
			mincInfo@volumeDataType <- dataType.df$numCode[enumCode]

			mincInfo@spaceType <- "talairach_"
			mincInfo@nDimensions <- 3
			# dim info ordering is zyx
			mincInfo@dimInfo <- data.frame(sizes=c(156, 220, 172), 
											steps=c(1, 1, 1), 
											starts=c(-68.25, -126.51, -86.095), 
											row.names=c("zspace", "yspace", "xspace"), 
											units=c("mm", "mm", "mm"), stringsAsFactors=FALSE)
		}


		# (3) volume should look like the MNI 305 volume used for PET
		if ( likeTemplate == "mni305PET" ) {

			# set the MincInfo fields
			mincInfo@filename <- filename

			dataClass.df <- rminc.getDataClasses()
			enumCode <- which(dataClass.df$string == "REAL")
			mincInfo@volumeDataClass <- dataClass.df$numCode[enumCode]

			dataType.df <- rminc.getDataTypes()
			enumCode <- which(dataType.df$code == "MI_TYPE_SHORT")
			mincInfo@volumeDataType <- dataType.df$numCode[enumCode]

			mincInfo@spaceType <- "talairach_"
			mincInfo@nDimensions <- 3
			# dim info ordering is zyx
			mincInfo@dimInfo <- data.frame(sizes=c(80, 128, 128), 
											steps=c(1.5, 1.72, 1.34), 
											starts=c(-37.5, -126.08, -85.76), 
											row.names=c("zspace", "yspace", "xspace"), 
											units=c("mm", "mm", "mm"), stringsAsFactors=FALSE)
		}


		# OK. We've now built a descriptive MincInfo object
		# new use this mincInfo to create a new, empty volume
		mincVolume <- mincIO.makeNewVolumeX(mincInfo)

		# DONE. Return MincVolumeIO object.
		if ( R_DEBUG_mincIO ) cat("<< mincIO.makeNewVolume()/templates ... \n")
		return(mincVolume)

		
	}
)



# Create a new, empty MincVolumeIO object, by passing a filename and a "like" volume name
setMethod(
	"mincIO.makeNewVolume", 
	signature=signature(filename="character", 
						likeTemplate="missing",
						likeFile="character"),
	definition=function(filename=filename, 
						dimLengths=NULL, dimSteps=NULL, dimStarts=NULL,
						likeTemplate=NULL,
						likeFile="dummy") {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.makeNewVolume()/like_file ... \n")

		# make sure that the "like" file exists, and is minc
		if ( !rminc.isMinc(likeFile) ) {
			stop(sprintf("Error: The \"like\" file [%s] either does not exist, or is not minc", likeFile));
		}
		
		# Ok, it's minc.  Make sure it's minc2
		likeFile <- rminc.asMinc2(likeFile);
		

		# Great, off to the races.  
		mincInfo <- mincIO.readMincInfo(likeFile);

		# update filename and then instantiate a MincVolumeIO object
		mincInfo@filename <- filename
		mincVolume <- mincIO.makeNewVolumeX(mincInfo)

		# DONE. Return MincVolumeIO object.
		if ( R_DEBUG_mincIO ) cat("<< mincIO.makeNewVolume()/like_file ... \n")
		return(mincVolume)

		
	}
)



mincIO.makeNewVolumeX <- function(mincInfo) {
#
# =============================================================================
# Purpose: Common 3D volume creation routine.
#
# Description:
#	Create and return a new MincVolumeIO object, given a MincInfo object. 
#
# TODO: 
#    (1) Nuttin.
#
#
# Notes: 
#
# =============================================================================
#
	if ( R_DEBUG_mincIO ) cat(">> mincIO.makeNewVolumeX() ... \n")

	# create a pre-initialized (to zero) 3D array
	nVoxels <- cumprod(mincInfo@dimInfo$sizes)[mincInfo@nDimensions]
	volumeData <- array(rep(0, nVoxels), dim=rev(mincInfo@dimInfo$sizes))

	# create the MincVolumeIO object
	mincVolume <- new("MincVolumeIO",
						volumeData,
						mincInfo=mincInfo, 
						volumeIntensityRange=mincInfo@volumeIntensityRange)

   # set display-related properties for an anatomical volume
	mincVolume@volumeType <- "anatomical"
	mincVolume@colorMap <- "gray"

	# DONE. Return MincVolumeIO object
	if ( R_DEBUG_mincIO ) cat("<< mincIO.makeNewVolumeX() ... \n")
	return(mincVolume)

}




# =============================================================================
# Purpose:	Given a 3D matrix, convert it to a MincVolumeIO object
#
#
# =============================================================================
# 
setGeneric( 
	name="mincIO.asVolume", 
	def = function(array3D,
					likeVolObject,
					likeTemplate) { standardGeneric("mincIO.asVolume") }
) 


# Convert to volume, by passing a MincVolumeIO object
setMethod(
	"mincIO.asVolume", 
	signature=signature(array3D="array", 
						likeVolObject="MincVolumeIO", 
						likeTemplate="missing"),
	definition=function(array3D, likeVolObject, likeTemplate=NULL) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.asVolume (likeVolume) ... \n")
		
		# Ok, so we're going to make the passed array into a volume
		# object that looks like the passed object
		
		# validation.  Make sure that the dimension sizes match.
		sameDim <- !all.equal(dim(likeVolObject), dim(array3D))
		if ( class(sameDim) == "character" ) {
			stop(sprintf("Array and \"like\" volume differ in dimensions. Conversion not permitted\n"))
		}

		# copy the "like" volume to get the S4 attributes
		myVolume <- likeVolObject

		# move the new array data in
		myVolume <- setDataPart(myVolume, array3D)
		
		# update volume-related fields
		myVolume@volumeIntensityRange <- c(min(myVolume), max(myVolume))
		myVolume@mincInfo@volumeIntensityRange <- c(myVolume@volumeIntensityRange[1],
			                                        myVolume@volumeIntensityRange[2])


		# DONE. Return MincInfo object.
		if ( R_DEBUG_mincIO ) cat("<< mincIO.asVolume (likeVolume) ... \n")
		return(myVolume)
		
	}
)


# Convert to volume, by passing a volume template type
setMethod(
	"mincIO.asVolume", 
	signature=signature(array3D="array", 
						likeVolObject="missing", 
						likeTemplate="character"),
	definition=function(array3D, likeVolObject=NULL, likeTemplate) {

		if ( R_DEBUG_mincIO ) cat(">> mincIO.asVolume (likeTemplate) ... \n")
		
		# Ok, so we're going to make the passed array into a volume
		# object that looks like the passed template volume
		
		# create the "like" volume to get the S4 attributes
		myVolume <- mincIO.makeNewVolume(filename="created_by_asVolume_function", likeTemplate=likeTemplate)

		# move the new array data in
		myVolume <- setDataPart(myVolume, array3D)
		
		# update volume-related fields
		myVolume@volumeIntensityRange <- c(min(myVolume), max(myVolume))
		myVolume@mincInfo@volumeIntensityRange <- c(myVolume@volumeIntensityRange[1],
			                                        myVolume@volumeIntensityRange[2])


		# DONE. Return MincInfo object.
		if ( R_DEBUG_mincIO ) cat("<< mincIO.asVolume (likeTemplate) ... \n")
		return(myVolume)
		
	}
)




