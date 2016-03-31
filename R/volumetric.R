

#' Explode a Label Volume into its Components
#' 
#' Given a label volume, this function splits the volume by label and then
#' returns a list() containing a mask volume for each of the labels.
#' 
#' @param label_vol A string containing the fully-qualified path to the input
#' label volume.
#' @param labels Options vector of label names
#' @param civetLabels A logical variable indicating whether the label volume is
#' using the Civet convention with regards to naming tissue type (e.g.,
#' 0=background, 1=csf, etc). If TRUE, the returned list components are named
#' using Civet tissue types (bg, csf, gm. wm), else components are simply
#' labelled by label number e.g. (``label_0'', ``label_2'', etc.).
#' @return A list is returned with each list item holding a mask volume
#' reflecting a particular label.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
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
		
		# compute the mask
		#cat(sprintf("processing label %s\n", label_name))
		mask_vol <- ifelse(label_vol == label, 1, 0)
		
		# adjust volume attributes to type of "mask"
		# mincIO.setProperty(mask_vol, "volumeType", "mask")
		# mincIO.setProperty(mask_vol, "colorMap", "gray")
		# mincIO.setProperty(mask_vol, "volumeIntensityRange", c(0,1))

		# save it to the return list
		return_list[[label_name]] <- mask_vol
	}
	
	return(return_list)
}


#' Combine Multiple Mask Volumes into a Single Mask Volume 
#' 
#' Given a list containing more than one label volume, combine those 
#' volumes, creating an aggregate mask volume. 
#' 
#' @param vol_list A list containing more than one mask volume.  Note that all 
#' volumes must reflect the same sampling. 
#' @return A single aggregate MincVolumeIO volume is returned. 
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
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
	#vol_sum <- mincIO.asVolume(vol_sum, vol_list[[1]])
	return(vol_sum)
}



