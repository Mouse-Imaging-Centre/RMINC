# =============================================================================
#
#     All S4 generic function definitions are here
#
# =============================================================================



# common S3 generics converted to S4
setGeneric("print") 
setGeneric("plot") 




# volume granularity generics
setGeneric( 
	name="mincIO.readVolume", 
	def = function(object, frameNo=0, ..., volumeType, colorMap) { standardGeneric("mincIO.readVolume") }
) 

setGeneric( 
	name="mincIO.writeVolume", 
	def = function(object, filename=missing, clobber=FALSE) { standardGeneric("mincIO.writeVolume") }
) 

setGeneric( 
	name="mincIO.makeNewVolume", 
	def = function(filename, 
					dimLengths, 
					dimSteps, 
					dimStarts,
					likeTemplate,
					likeFile) { standardGeneric("mincIO.makeNewVolume") }
) 

setGeneric( 
	name="mincIO.asVolume", 
	def = function(array3D,
					likeVolObject,
					likeTemplate) { standardGeneric("mincIO.asVolume") }
) 



# slice granularity generics
setGeneric( 
	name="mincIO.getSliceX", 
	def = function(mincVolume, sliceNo) { standardGeneric("mincIO.getSliceX") }
) 

setGeneric( 
	name="mincIO.getSliceY", 
	def = function(mincVolume, sliceNo) { standardGeneric("mincIO.getSliceY") }
) 

setGeneric( 
	name="mincIO.getSliceZ", 
	def = function(mincVolume, sliceNo) { standardGeneric("mincIO.getSliceZ") }
) 

setGeneric( 
	name="mincIO.putSlice", 
	def = function(mincSlice, mincVolume, sliceNo, ...) { standardGeneric("mincIO.putSlice") }
) 

setGeneric( 
	name="mincIO.makeNewSliceZ", 
	def = function(mincVolume, initVector=missing) { standardGeneric("mincIO.makeNewSliceZ") }
) 


# sliceIO generics
setGeneric( 
	name="mincIO.readBySlice", 
	def = function(filenames, sliceNumber, ..., volumeType, colorMap) { standardGeneric("mincIO.readBySlice") }
) 

setGeneric( 
	name="mincIO.getSliceFromSliceArray", 
	def = function(mincSliceMatrix, sliceIndex) { standardGeneric("mincIO.getSliceFromSliceArray") }
) 


# 
# voxelIO generics
setGeneric( 
	name="mincIO.readByVoxel", 
	def = function(filenames, voxelCoords) { standardGeneric("mincIO.readByVoxel") }
) 








