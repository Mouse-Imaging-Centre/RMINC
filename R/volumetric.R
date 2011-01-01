

# =============================================================================
# Purpose: 
#	Given a labelled MincVolumeIO object, split volume into separate volumes,
#	one for each label type.  Return as a list of MincVolumeIO objects.
#
# Example:
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	labelVolname <- civet.getFilenameClassify(scanID, baseDir)
#	label_vol <- mincIO.readVolume(labelVolname)
#	components <- volume.explodeLabelVolume(label_vol)
#
# Note:
#	This function returns a list of mincIO objects, with each object
#	labelled either by its label number, or by the tissue type name
#	(if civetLabels == TRUE).
# =============================================================================
#
volume.explodeLabelVolume <- function(label_vol, labels=NULL, civetLabels=TRUE) {
	
	# list the tissue types in Civet order, to be used in naming the
	# components of the returned list
	tissueTypes <- c("bg", "csf", "gm", "wm")
	
	# read the label volume
	# label_vol <- mincIO.readVolume(labelVolname)

	# here we have a choice: We either specify which labels we want exploded out,
	# or, we don't, and we get all of them
	#
	if ( is.null(labels)) {
		# not specified, so get the set of unique labels
		labels <- unique(label_vol)
	} else {
		# specified, ... make sure we've been passed a numeric vector
		if ( !is.vector(labels, mode="numeric")) {stop("Label vector must be a numeric vector")}
		
	}


	# loop thru all labels, creating a mask for each
	return_list <- list()
	for ( label in labels ) {
		#
		# determine the list component name
		label_name <- paste("label", label, sep="_")
		if ( civetLabels ) { label_name <- tissueTypes[label +1]}
		if (R_DEBUG_mincIO) cat(sprintf("processing label %s\n", label_name))
		
		# compute the mask
		#cat(sprintf("processing label %s\n", label_name))
		mask_vol <- ifelse(label_vol == label, 1, 0)
		
		# adjust volume attributes to type of "mask"
		mincIO.setProperty(mask_vol, "volumeType", "mask")
		mincIO.setProperty(mask_vol, "colorMap", "gray")
		mincIO.setProperty(mask_vol, "volumeIntensityRange", c(0,1))

		# save it to the return list
		return_list[[label_name]] <- mask_vol
	}
	
	return(return_list)
}


# =============================================================================
# Purpose: 
#	Given a list binarized MincVolumeIO objects, combine them into one,
#	returning the single, combined MincVolumeIO mask object.
#
# Example:
#	labelVolname <- civet.getFilenameClassify(scanID, baseDir)
#	label_vol <- mincIO.readVolume(labelVolname)
#	vol_list <- volume.explodeLabelVolume(label_vol)
#
#	# ... remove background component -- combine the rest
#	xl <- list(vol_list$gm, vol_list$wm, vol_list$csf)
#	cVol <- volume.combineMaskVolumes(xl)
#
# Note: Ummmm ... nothing much.
# =============================================================================
#
volume.combineMaskVolumes <- function(vol_list) {

	# first, get the number of volumes to combine
	nVolumes <- length(vol_list)

	# loop thru all labels, adding 'em up
	vol_sum <- vol_list[[1]]
	for ( vol in 2:nVolumes ) {
		vol_sum <- vol_sum + vol_list[[vol]]
	}
	
	# set voxels included in multiple masks to '1'
	vol_sum <- ifelse(round(vol_sum) < 1, 0, 1)
	#
	# make sure it's still a volume, then return
	vol_sum <- mincIO.asVolume(vol_sum, vol_list[[1]])
	return(vol_sum)
}



