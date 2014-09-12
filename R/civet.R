
# =============================================================================
# Purpose: 
#	Check whether the passed Civet version is one for which we have tested.
#
# Example:
#	civet.checkVersion("1.1.8")
# =============================================================================
#
civet.checkVersion <- function(civetVersion) {
	if ( civetVersion != "1.1.9" &&  civetVersion != "1.1.7") {
		warning(sprintf("This function has not been tested with Civet version %s. Use at your own risk.", civetVersion), immediate.=TRUE)
	}
	return
}


# =============================================================================
# Purpose: 
#	Return the fully-qualified filename of the *_classify.mnc volume
#
# Example:
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	filename <- civet.getFilenameClassify(scanID, baseDir)
#
# Note:
#	This class of functions actually finds the filename by looking for
#	it in the appropriate Civet subdir. As such, the file must actually
#	exist for this function to work.
# =============================================================================
#
civet.getFilenameClassify <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# create the scan-level root dir name
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	
	# get a list of matching filenames in the classify dir, and return
	classifyDir <- file.path(scanRoot, 'classify')
	files <- list.files(classifyDir, pattern=glob2rx("*_classify.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(classifyDir, files)
	}
	
	return(files)
}


# =============================================================================
# Purpose: 
#	Return the fully-qualified filename of the *_pve_*.mnc volumes
#
# Example: See civet.getFilenameClassify
# Note: civet.getFilenameClassify
# =============================================================================
#
civet.getFilenameGrayMatterPve <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	classifyDir <- file.path(scanRoot, 'classify')
	files <- list.files(classifyDir, pattern=glob2rx("*_pve_gm.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(classifyDir, files)
	}
	
	return(files)
}
	#
civet.getFilenameWhiteMatterPve <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	classifyDir <- file.path(scanRoot, 'classify')
	files <- list.files(classifyDir, pattern=glob2rx("*_pve_wm.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(classifyDir, files)
	}
	
	return(files)
}
	#
civet.getFilenameCsfPve <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	classifyDir <- file.path(scanRoot, 'classify')
	files <- list.files(classifyDir, pattern=glob2rx("*_pve_csf.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(classifyDir, files)
	}
	
	return(files)

}


# =============================================================================
# Purpose: 
#	Return the fully-qualified filename of the *_t1_final.mnc volumes
#
# Example: See civet.getFilenameClassify
# Note: civet.getFilenameClassify
# =============================================================================
#
civet.getFilenameStxT1 <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	finalDir <- file.path(scanRoot, 'final')
	files <- list.files(finalDir, pattern=glob2rx("*_t1_final.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(finalDir, files)
	}
	
	return(files)
}


# =============================================================================
# Purpose: 
#	Return the fully-qualified filename of the *_brain_mask.mnc or 
#	*_skull_mask.mnc volumes.
#
# Example: See civet.getFilenameClassify
# Note: The brain mask differentiates itself from the skull mask in that
#		the brain mask does *not* include the cerebellum.
# =============================================================================
#
civet.getFilenameCerebrumMask <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	maskDir <- file.path(scanRoot, 'mask')
	files <- list.files(maskDir, pattern=glob2rx("*_brain_mask.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(maskDir, files)
	}
	
	return(files)
}
#
civet.getFilenameSkullMask <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	maskDir <- file.path(scanRoot, 'mask')
	files <- list.files(maskDir, pattern=glob2rx("*_skull_mask.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(maskDir, files)
	}
	
	return(files)
}


# =============================================================================
# Purpose: 
#	Return the fully-qualified gray/mid/white matter surface filenames 
#	(left and right).
#
# Example: See civet.getFilenameClassify
# Note: civet.getFilenameClassify
# =============================================================================
#
civet.getFilenameGrayMatterSurfaces <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	surfacesDir <- file.path(scanRoot, 'surfaces')
	
	# do we have rsl files? (we can get this in earlier Civet versions)
	files <- list.files(surfacesDir, pattern=glob2rx("*_gray_surface_rsl_*.obj"))
	# ... if not, look for non-resampled surfaces 
	if ( length(files) == 0 ) {
		files <- list.files(surfacesDir, pattern=glob2rx("*_gray_surface_*.obj"))
	}
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(surfacesDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}
#
civet.getFilenameWhiteMatterSurfaces <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	surfacesDir <- file.path(scanRoot, 'surfaces')
	
	# do we have rsl files? (we can get this in earlier Civet versions)
	files <- list.files(surfacesDir, pattern=glob2rx("*_white_surface_rsl_*_calibrated_*.obj"))
	# ... if not, look for non-resampled surfaces 
	if ( length(files) == 0 ) {
		files <- list.files(surfacesDir, pattern=glob2rx("*_white_surface_*_calibrated_*.obj"))
	}
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(surfacesDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}
#
civet.getFilenameMidSurfaces <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	surfacesDir <- file.path(scanRoot, 'surfaces')
	
	# do we have rsl files? (we can get this in earlier Civet versions)
	files <- list.files(surfacesDir, pattern=glob2rx("*_mid_surface_rsl_*.obj"))
	# ... if not, look for non-resampled surfaces 
	if ( length(files) == 0 ) {
		files <- list.files(surfacesDir, pattern=glob2rx("*_mid_surface_*.obj"))
	}
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(surfacesDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}


# =============================================================================
# Purpose: 
#	Return the fully-qualified cortical thickness filenames 
#	(left and right).
#
# Example: See civet.getFilenameClassify
# Note: civet.getFilenameClassify
# =============================================================================
#
civet.getFilenameCorticalThickness <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	ctDir <- file.path(scanRoot, 'thickness')
	files <- list.files(ctDir, pattern=glob2rx("*_native_rms_rsl_tlink_20mm_*.txt"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(ctDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}
#

# =============================================================================
# Purpose: 
#	Return the fully-qualified mean surface curvature filenames 
#	(left and right).
#
# Example: See civet.getFilenameClassify
# Note: civet.getFilenameClassify
# =============================================================================
#
civet.getFilenameMeanSurfaceCurvature <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	ctDir <- file.path(scanRoot, 'thickness')
	files <- list.files(ctDir, pattern=glob2rx("*_native_mc_rsl_20mm_*.txt"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(ctDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}
#


# =============================================================================
# Purpose: 
#	Return the fully-qualified transform filenames (linear and nonlinear).
#
# Example: See civet.getFilenameClassify
# Note: The request for the nonlinear transform file, returns 2 filenames: 
#		the xfm filename and the grid volume name.
# =============================================================================
#
civet.getFilenameLinearTransform <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	xfmDir <- file.path(scanRoot, 'transforms/linear')
	files <- list.files(xfmDir, pattern=glob2rx("*_t1_tal.xfm"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(xfmDir, files)
	}
	
	return(files)
}
#
#
civet.getFilenameNonlinearTransform <- function(scanID, baseDir, civetVersion="1.1.9", fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	xfmDir <- file.path(scanRoot, 'transforms/nonlinear')
	file.xfm <- list.files(xfmDir, pattern=glob2rx("*_nlfit_It.xfm"))
	file.grid <- list.files(xfmDir, pattern=glob2rx("*_nlfit_It_grid_*.mnc"))
	
	# fully-qualified path requested?
	if ( fullPath ) {
		file.xfm <- file.path(xfmDir, file.xfm)
		file.grid <- file.path(xfmDir, file.grid)
	}

	return(list(xfm=file.xfm, grid=file.grid))
}


###########################################################################################
#' @description Generates list of filenames output by CIVET
#' @name civet.getAllFilenames
#' @title civet.getAllFilenames
#' @usage civet.getAllFilenames(gf, idvar, prefix, basedir, append = TRUE, civetVersion = "1.1.9")
#' @param gf Data Frame with subject information
#' @param idvar column name in gf with subject IDs
#' @param prefix Prefix specified when CIVET was run
#' @param basedir directory where all CIVET output was stored
#' @param append Whether to append the results to the input gf
#' @param civetVersion Version of CIVET 
#' @details Prior to running, read.csv  may be called to generate the input argument gf. 
#' The results will be stored under the column name CIVETFILES either in the input gf (if append = TRUE) or in a new gf. 
#' Currently only CIVET versions 1.1.9 and 1.1.12 are supported.
#' @return gf is returned with CIVET filenames 
#' @seealso civet.readAllCivetFiles
#' @examples
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
###########################################################################################

civet.getAllFilenames <- function(gf, idvar, prefix, basedir, append=TRUE, civetVersion="1.1.9") {
	# designed for use with CIVET 1.1.9 and CIVET 1.1.12
	if ( civetVersion != "1.1.9"  && civetVersion != "1.1.12" ) {
		warning("This function has only been tested with Civet version 1.1.9. and 1.1.12 Use at your own risk.", immediate.=TRUE)
	}
	
	# extract the scanIDs from the glim
	ids <- gf[,idvar]
	
	# create the scan-level root dir names
	b <- paste(basedir, "/", ids, "/", sep="")

	# insert fully-qualified file names
	filenames.df <- data.frame(leftGIFiles=rep(0,nrow(gf)))

	if (civetVersion == "1.1.12")
	{	
  		filenames.df$leftGIFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_gi_left.dat",sep = "")
  		filenames.df$rightGIFiles <- paste(b, "surfaces/", prefix, "_", ids, "_gi_right.dat",sep = "")
  		filenames.df$leftlobeArea40mmFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_areas_40mm_left.dat",sep = "")
  		filenames.df$rightlobeArea40mmFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_areas_40mm_right.dat",sep = "")
  		filenames.df$leftlobeThicknessFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_thickness_tlink_20mm_left.dat",sep = "")
  		filenames.df$rightlobeThicknessFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_thickness_tlink_20mm_right.dat",sep = "")
  		filenames.df$leftlobeVolumeFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_volumes_40mm_left.dat",sep = "")
  		filenames.df$rightlobeVolumeFiles  <- paste(b, "surfaces/", prefix, "_", ids, "_lobe_volumes_40mm_right.dat",sep = "")
		filenames.df$midSurfaceleftNativeArea = paste(b, "surfaces/", prefix, "_", ids, "_mid_surface_rsl_left_native_area_40mm.txt",sep = "")
		filenames.df$midSurfacerightNativeArea = paste(b,"surfaces/", prefix, "_", ids, "_mid_surface_rsl_right_native_area_40mm.txt",sep = "")
		filenames.df$SurfaceleftNativeVolume = paste(b, "surfaces/", prefix, "_", ids, "_surface_rsl_left_native_volume_40mm.txt",sep = "")
		filenames.df$SurfacerightNativeVolume = paste(b,"surfaces/", prefix, "_", ids, "_surface_rsl_right_native_volume_40mm.txt",sep = "") 
 
  		filenames.df$brain_volume = paste(b, "classify/", prefix, "_", ids, "_cls_volumes.dat",sep = "")
		filenames.df$cerebral_volume = paste(b, "thickness/", prefix, "_", ids, "_cerebral_volume.dat",sep = "")
		filenames.df$nativeRMS_RSLtlink20mmleft = paste(b,"thickness/", prefix, "_", ids, "_native_rms_rsl_tlink_20mm_left.txt",sep = "") 
		filenames.df$nativeRMS_RSLtlink20mmright = paste(b,"thickness/", prefix, "_", ids, "_native_rms_rsl_tlink_20mm_right.txt",sep = "") 
		filenames.df$nativeRMStlink20mmleft = paste(b,"thickness/", prefix, "_", ids, "_native_rms_tlink_20mm_left.txt",sep = "") 
		filenames.df$nativeRMStlink20mmright = paste(b,"thickness/", prefix, "_", ids, "_native_rms_tlink_20mm_right.txt",sep = "") 



		gf$CIVETFILES = filenames.df;
		return(gf)

	}
	
	else

	{
		filenames.df$tissue <- paste(b, "classify/", prefix, "_", ids, "_cls_volumes.dat",
		           sep="")
		filenames.df$structures <- paste(b, "segment/", prefix, "_", ids, "_masked.dat",
		               sep="")
		filenames.df$left.thickness <- paste(b, "thickness/", prefix, "_", ids,
		                   "_native_rms_rsl_tlink_20mm_left.txt", sep="")
		filenames.df$right.thickness <- paste(b, "thickness/", prefix, "_", ids,
		                   "_native_rms_rsl_tlink_20mm_right.txt", sep="")
		filenames.df$rightROIthickness <- paste(b, "segment/", prefix, "_", ids,
		                      "_lobe_thickness_tlink_20mm_right.dat", sep="")
		filenames.df$leftROIthickness <- paste(b, "segment/", prefix, "_", ids,
		                     "_lobe_thickness_tlink_20mm_left.dat", sep="")
		filenames.df$rightROIarea <- paste(b, "segment/", prefix, "_", ids,
		                 "_lobe_areas_right.dat", sep="")
		filenames.df$leftROIarea <- paste(b, "segment/", prefix, "_", ids,
		                "_lobe_areas_left.dat", sep="")
		filenames.df$GMVBM <- paste(b, "VBM/", prefix, "_", ids,
		          "_smooth_8mm_gm.mnc", sep="")
		filenames.df$WMVBM <- paste(b, "VBM/", prefix, "_", ids,
		          "_smooth_8mm_wm.mnc", sep="")
		filenames.df$CSFVBM <- paste(b, "VBM/", prefix, "_", ids,
		           "_smooth_8mm_csf.mnc", sep="")
		filenames.df$GMVBMsym <- paste(b, "VBM/", prefix, "_", ids,
		          "_smooth_8mm_gm_sym.mnc", sep="")
		filenames.df$WMVBMsym <- paste(b, "VBM/", prefix, "_", ids,
		          "_smooth_8mm_wm_sym.mnc", sep="")
		filenames.df$CSFVBMsym <- paste(b, "VBM/", prefix, "_", ids,
		             "_smooth_8mm_csf_sym.mnc", sep="")
  	}
	# append names to the input glim, if desired
	if ( append == TRUE) {
		filenames.df <- cbind(gf, filenames.df)
	}

	return(filenames.df)
}




civet.AllROIs <- function(gf, defprefix) {
	# =============================================================================
	# Purpose: Ummmmm ..... 
	#
	# Note: Original code written by Jason Lerch, I (Jim) just changed the name
	#		to use the new civet.* prefix.
	#
	# =============================================================================
	#
  tissues <- anatGetAll(gf$tissue, NULL, method="text",
                        defs=paste(defprefix,"/","tissue-volumes.csv",sep=""),
                        drop=TRUE)
  colnames(tissues) <- paste(colnames(tissues), "volume")

  structures <- anatGetAll(gf$structures, NULL, method="text",
                           defs=paste(defprefix, "/","lobe-seg-volumes.csv",
                             sep=""),
                           drop=TRUE)
  colnames(structures) <- paste(colnames(structures), "volume")
  leftareas <- anatGetAll(gf$leftROIarea, NULL, method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="left")
  colnames(leftareas) <- paste(colnames(leftareas), "surface area")
  rightareas <- anatGetAll(gf$rightROIarea, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="right")
  colnames(rightareas) <- paste(colnames(rightareas), "surface area")
  leftthickness <- anatGetAll(gf$leftROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="left")
  colnames(leftthickness) <- paste(colnames(leftthickness), "mean thickness")
  rightthickness <- anatGetAll(gf$rightROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="right")
  colnames(rightthickness) <- paste(colnames(rightthickness), "mean thickness")

  vols <- data.frame(tissues, structures, leftareas, rightareas, leftthickness, rightthickness)
  return(vols)
}
# =============================================================================
# Purpose: 
#	Organizes CIVET .dat file based on atlas
#
# Example:
#	atlasFile <- "AAL.csv"
#	dataFiles <- list of data files to organize
#
# =============================================================================
civet.organizeCivetDatFilesAtlas <- function(atlasFile,dataFiles, civetVersion="1.1.12") {
	
	#Initializaion
	numberOfFiles = length(dataFiles)
	gf = read.csv(atlasFile)
	roiTable = matrix(data=NA,nrow=length(gf[,1]),ncol=numberOfFiles/2)
	roiLabels = c()

	# Label rows with ROI labels
	for (i in 1:length(gf[,1])) { roiLabels[i] = as.character(gf[i,3])}
	rownames(roiTable) = roiLabels


	# Iterate through each file, look up ROI location and fill in  table
	for (j in 1:numberOfFiles)
	{
		halfFiles = numberOfFiles / 2;
		if(j > halfFiles)
		{
			columnIndex = j - numberOfFiles/2
		}
		else
		{	
			columnIndex = j;
		}
		if(file.exists(dataFiles[j]))
		{
			labels = read.table(dataFiles[j])[,1]
			value = read.table(dataFiles[j])[,2]
			labels = as.numeric(as.character(labels))
			value = as.numeric(as.character(value))
			for (i in 1:length(labels)) 
			{ 
				roiName = as.character( gf[which(gf == labels[i],arr.ind=FALSE),3])
				tableIndex = which(rownames(roiTable) == roiName)
				roiTable[tableIndex,columnIndex] = value[i]
			}
		}	
	}

	return(roiTable)
}
# =============================================================================
# Purpose: 
#  Organizes CIVET .dat file based on whole brain,grey and white volumes
#
# Example:
#	dataFiles <- list of data files to organize
#
# =============================================================================
civet.organizeCivetDatFilesWholeBrain<- function(dataFiles, civetVersion="1.1.12") {
  
  numberOfFiles = length(dataFiles)
  roiTable = matrix(data=NA,nrow=3,ncol=numberOfFiles)
  rownames(roiTable) = c("CSF","GM","WM")
  for (j in 1:numberOfFiles)
  {
    
    if(file.exists(dataFiles[j]))
    {
      value = read.table(dataFiles[j])[,2]
      value = as.numeric(as.character(value))
      roiTable[1:3,j] = value
    }	
    
  }
  return(roiTable)
}

# =============================================================================
# Purpose: 
#	Organizes CIVET .dat file based on mid, white and grey labels
#
# Example:
#	dataFiles <- list of data files to organize
# =============================================================================
civet.organizeCivetDatFilesMidWhiteGrey <- function(dataFiles, civetVersion="1.1.12") {

	numberOfFiles = length(dataFiles)
	roiTable = matrix(data=NA,nrow=6,ncol=numberOfFiles/2)
	rownames(roiTable) = c("leftGray", "leftWhite","leftMid","rightGray","rightWhite","rightMid")

	for (j in 1:numberOfFiles)
	{	
		halfFiles = numberOfFiles / 2;
		if(j > halfFiles)
		{
			columnIndex = j - numberOfFiles/2
			rowIndex = 4;
		}
		else
		{	
			columnIndex = j;
			rowIndex = 1;
		}

		if(file.exists(dataFiles[j]))
		{
			value = read.table(dataFiles[j])[,4]
			value = as.numeric(as.character(value))
			if(length(value) == 3)
				roiTable[rowIndex:(rowIndex+2),columnIndex] = value
		}	

	}
	return(roiTable)
}
# =============================================================================
# Purpose: 
#	Organizes CIVET .txt file with vertex measures
#
# Example:
#	dataFiles <- list of data files to organize
#
# Note : Left Vertices are listed First
# =============================================================================
civet.organizeCivetTxtFilesVertex <- function(dataFiles) {

	numberOfFiles = length(dataFiles)
	roiTable = matrix(data=NA,nrow=40962*2,ncol=numberOfFiles/2)

	for (j in 1:numberOfFiles)
	{
		halfFiles = numberOfFiles / 2;
		if(j > halfFiles)
		{
			columnIndex = j - numberOfFiles/2
			rowIndex = 40963;
		}
		else
		{	
			columnIndex = j;
			rowIndex = 1;
		}

		if(file.exists(dataFiles[j]))
		{
			value = read.table(dataFiles[j])[,1]
			value = as.numeric(as.character(value))
			roiTable[rowIndex:(rowIndex+40961),columnIndex] = value
		}	

	}
	return(roiTable)
}

###########################################################################################
#' @description Parses outputs from CIVET pipeline 
#' @name civet.readAllCivetFiles
#' @title Read  all CIVET files into R
#' @usage civet.readAllCivetFiles(atlasFile, gf)
#' @param gf Data Frame containing list of all CIVET file names, and where results will be stored
#' @param atlasFile Full path to the atlas on which CIVET was run ;
#' @details Prior to running, civet.getAllFilenames may be called to generate the input argument gf .
#' This function will extract the following information from the CIVET pipeline: Lobe Area (40 mm),
#' Lobe Thickness, Lobe Volume, GI, Mid Surface Native Area, Surface Native Volume, Native RMS RSL tLink (20mm), Native RMS tLink (20 mm)
#' @return gf is returned with CIVET values
#' @seealso  civet.getAllFilenames
#' @examples
#' getRMINCTestData() 
#' gf = read.csv("~/SubjectTable.csv") 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
###########################################################################################
civet.readAllCivetFiles = function(atlasFile,gf)
{
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Lobe Area, Lobe Thickness, Lobe Volume
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	roiTable = civet.organizeCivetDatFilesAtlas(atlasFile,c(gf$CIVETFILES$leftlobeArea40mmFiles, gf$CIVETFILES$rightlobeArea40mmFiles))
	gf$lobeArea40mm = t(roiTable)

	roiTable = civet.organizeCivetDatFilesAtlas(atlasFile,c(gf$CIVETFILES$leftlobeThicknessFiles, gf$CIVETFILES$rightlobeThicknessFiles))
	gf$lobeThickness = t(roiTable)

	roiTable = civet.organizeCivetDatFilesAtlas(atlasFile,c(gf$CIVETFILES$leftlobeVolumeFiles, gf$CIVETFILES$rightlobeVolumeFiles))
	gf$lobeVolume = t(roiTable)

	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# GI
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	roiTable = civet.organizeCivetDatFilesMidWhiteGrey(c(gf$CIVETFILES$leftGIFiles, gf$CIVETFILES$rightGIFiles))	
	gf$GI = t(roiTable)
	
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Cerebral Volume
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
	#roiTable = civet.organizeCivetDatFilesWholeBrain(c(gf$CIVETFILES$cerebral_volume))  
	#gf$cerebralVolume = t(roiTable)
  
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =civet.readAllCivetFiles
	# Brain Volume
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
	roiTable = civet.organizeCivetDatFilesWholeBrain(c(gf$CIVETFILES$brain_volume))  
	gf$BrainVolume = t(roiTable)
  
  
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Vertex based Measures
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	roiTable = civet.organizeCivetTxtFilesVertex(c(gf$CIVETFILES$midSurfaceleftNativeArea, gf$CIVETFILES$midSurfacerightNativeArea))
	gf$midSurfaceNativeArea = t(roiTable)
	
	roiTable = civet.organizeCivetTxtFilesVertex(c(gf$CIVETFILES$SurfaceleftNativeVolume, gf$CIVETFILES$SurfacerightNativeVolume))
	gf$SurfaceNativeVolume = t(roiTable)

	roiTable = civet.organizeCivetTxtFilesVertex(c(gf$CIVETFILES$nativeRMS_RSLtlink20mmleft, gf$CIVETFILES$nativeRMS_RSLtlink20mmright))
	gf$nativeRMS_RSLtlink20mm = t(roiTable)

	roiTable = civet.organizeCivetTxtFilesVertex(c(gf$CIVETFILES$nativeRMStlink20mmleft, gf$CIVETFILES$nativeRMStlink20mmright))
	gf$nativeRMStlink20mm= t(roiTable)

	return(gf)
}


# =============================================================================
# Purpose: 
#	Read a selection of Civet-generated *.dat files and return the
#	contents in a list.
#
# Example:
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	civet_dat_list <- civet.readCivetDatFiles(scanID, baseDir)
#
# Note: Nothing really.
# =============================================================================
#
civet.readCivetDatFiles <- function(scanID, baseDir, civetVersion="1.1.9") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# create the scan-level root dir name
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	
	# init returning list
	return_list <-list()
	
	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# load the tissue classify volume info
	#	
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'classify')
	myFile <- list.files(subDir, pattern=glob2rx("*_cls_volumes.dat"))
	myFile_fullPath <- file.path(subDir, myFile)
	myDf <- read.table(myFile_fullPath,
	row.names=c("csf", "gm", "wm"),
	col.names=c("tissueType_code", "volume"))
	#	
	return_list$native_tissue_volumes <- myDf



	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# load gyrification indices
		
	
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'surfaces')
	# ... left
	gi_left_file <- list.files(subDir, pattern=glob2rx("*_gi_left.dat"))
	gi_left_file_fullPath <- file.path(subDir, gi_left_file)
	gi_left.df <- read.table(gi_left_file_fullPath)
	# ... right
	gi_right_file <- list.files(subDir, pattern=glob2rx("*_gi_right.dat"))
	gi_right_file_fullPath <- file.path(subDir, gi_right_file)
	gi_right.df <- read.table(gi_right_file_fullPath)
	#
	return_list$gyrification_index <- c(lh=gi_left.df[1,3], rh=gi_right.df[1,3]) 



	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# load the cerebral volume info
	#	
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'thickness')
	cerebral_volumes <- list.files(subDir, pattern=glob2rx("*_cerebral_volume.dat"))
	cerebral_volumes_fullPath <- file.path(subDir, cerebral_volumes)
	cerebral_volumes.df <- read.table(cerebral_volumes_fullPath, 
	row.names=c("extra_cerebral_csf", "cortical_gray", "wmSurface_plus_contents"),
	col.names=c("code", "volume"))
	#	
	return_list$native_cerebral_volumes <- cerebral_volumes.df



	# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# load the percentage values from the verify subdir
	#	
	#
	
	# create containing list
	quality_control <- list()
	
	# (A) Load *_brainmask_qc.txt
	#	e.g., "native skull mask in stx space (11.77%)"
	#
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'verify')
	verify_qc <- list.files(subDir, pattern=glob2rx("*_brainmask_qc.txt"))
	verify_qc_fullPath <- file.path(subDir, verify_qc)
	# ... open file connection, read 1 line, then close()
	verify_qc.con <- file(verify_qc_fullPath, open="rt")
	txtLine <- readLines(verify_qc.con, n=1)
	close(verify_qc.con)
	#
	# ... extract the percentage value
	ndx <- regexpr("[[:digit:]]+.[[:digit:]]+", txtLine)
	pctValue <- substr(txtLine, ndx, ndx-1+attr(ndx, "match.length"))
	pctValue <- as.numeric(pctValue)
	#	
	quality_control$brainmask_qc <- pctValue


	#
	# (B) Load *_classify_qc.txt
	#	e.g., "classified image CSF 15.01%  GM 47.86%  WM 37.13%"
	#
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'verify')
	verify_qc <- list.files(subDir, pattern=glob2rx("*_classify_qc.txt"))
	verify_qc_fullPath <- file.path(subDir, verify_qc)
	# ... open file connection, read 1 line, then close()
	verify_qc.con <- file(verify_qc_fullPath, open="rt")
	txtLine <- readLines(verify_qc.con, n=1)
	close(verify_qc.con)
	
	
	# init receiving data.frame
	classify.df <- data.frame(pct=rep(0,3), row.names=c("csf", "gm", "wm"))
	
	# ... extract the indices to the percentage values
	ndx <- gregexpr("[[:digit:]]+.[[:digit:]]+", txtLine)
	csf_ndx <- c(ndx[[1]][1], 
	             ndx[[1]][1] - 1 + attr(ndx[[1]], "match.length")[1])

	gm_ndx <- c(ndx[[1]][2], 
	             ndx[[1]][2] - 1 + attr(ndx[[1]], "match.length")[2])

	wm_ndx <- c(ndx[[1]][3], 
	             ndx[[1]][3] - 1 + attr(ndx[[1]], "match.length")[3])
		
	# ... extract and store the values
	csf_pct <- substr(txtLine, csf_ndx[1], csf_ndx[2])
	classify.df["csf", "pct"] <- as.numeric(csf_pct)
	#	
	gm_pct <- substr(txtLine, gm_ndx[1], gm_ndx[2])
	classify.df["gm", "pct"] <- as.numeric(gm_pct)
	#
	wm_pct <- substr(txtLine, wm_ndx[1], wm_ndx[2])
	classify.df["wm", "pct"] <- as.numeric(wm_pct)
	#
	quality_control$classify_qc <- classify.df



	#
	# (C) Load *_surface_qc.txt
	#	e.g., "white surface ( 8.06%), gray surface (10.03%)"
	#
	# get a list of matching filenames in the classify dir, and return
	subDir <- file.path(scanRoot, 'verify')
	verify_qc <- list.files(subDir, pattern=glob2rx("*_surface_qc.txt"))
	verify_qc_fullPath <- file.path(subDir, verify_qc)
	# ... open file connection, read 1 line, then close()
	verify_qc.con <- file(verify_qc_fullPath, open="rt")
	txtLine <- readLines(verify_qc.con, n=1)
	close(verify_qc.con)
	
	
	# init receiving data.frame
	surface.df <- data.frame(pct=rep(0,2), row.names=c("gm", "wm"))
	
	# ... extract the indices to the percentage values
	ndx <- gregexpr("[[:digit:]]+.[[:digit:]]+", txtLine)
	wm_ndx <- c(ndx[[1]][1], 
	             ndx[[1]][1] - 1 + attr(ndx[[1]], "match.length")[1])

	gm_ndx <- c(ndx[[1]][2], 
	             ndx[[1]][2] - 1 + attr(ndx[[1]], "match.length")[2])
		
	# ... extract and store the values
	gm_pct <- substr(txtLine, gm_ndx[1], gm_ndx[2])
	surface.df["gm", "pct"] <- as.numeric(gm_pct)
	#	
	wm_pct <- substr(txtLine, wm_ndx[1], wm_ndx[2])
	surface.df["wm", "pct"] <- as.numeric(wm_pct)
	#
	quality_control$surface_qc <- surface.df
	
	#
	return_list$quality_control <- quality_control

	# return final list
	return(return_list)

}



# =============================================================================
# Purpose: 
#	Return a named vector containing the tisse volumes, in stx space, derived from 
#	the final tissue classification volume.
#
# Example: 
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	cls_vec <- civet.computeStxTissueVolumes(scanID, baseDir)
#	print(cls_vec["gm"])
# =============================================================================
#
civet.computeStxTissueVolumes <- function(scanID, baseDir, civetVersion="1.1.9") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)

	# get classify volume name
	filename <- civet.getFilenameClassify(scanID, baseDir)
	cls_vol <- mincIO.readVolume(filename)

	# explode classify into components
	clsX <- volume.explodeLabelVolume(cls_vol, civetLabels=TRUE)
	
	# store elements into named vector
	cls_vec <- numeric(3)
	names(cls_vec) <- c("csf", "gm", "wm")
	cls_vec["csf"] <- sum(clsX$csf)
	cls_vec["gm"] <- sum(clsX$gm)
	cls_vec["wm"] <- sum(clsX$wm)
	
	return(cls_vec)
}
#

# =============================================================================
# Purpose: 
#	XFM files contain information to transform a given volume to a model.
#	In the case of Civet and rescaling, the XFM contains the rescaling 
#	factors (x,y,z) needed to transform the Native volume to the model, which
#	currently, is usually the symmetrical icbm-152 model.
#
#	This functuon serves to compute a global rescaling factor by reading
#	the individual x,y,z rescales from the XFM, and returning the
#	product.
#
#	Interpretation of rescaling factors:
#	(a) > 1.0 = native brain is expanded to fit model
#	(b) < 1.0 = native brain is reduced to fit model
#
# Example: 
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	rescale <- civet.computeNativeToStxRescalingFactor(scanID, baseDir)
#	print(rescale)
# =============================================================================
#
civet.computeNativeToStxRescalingFactor <- function(scanID, baseDir, civetVersion="1.1.9") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# read linear xfm file
	xfmFilename <- civet.getFilenameLinearTransform(scanID, baseDir)
	xfm.df <- rminc.readLinearXfmFile(xfmFilename)
	
	# compute global linear scaling factor
	scaling_factor <- prod(xfm.df["scale",])
	
	# return
	return(scaling_factor)
}
#



# =============================================================================
# Purpose: 
#	Return a named vector containing NATIVE space tisse volumes, derived from 
#	the final tissue classification volume.
#
# Example: 
#	scanID <- "0548-F-NC"
#	baseDir <- "~/tmp/ADNI/civet/pipeOut"
#	cls_vec <- civet.computeNativeTissueVolumes(scanID, baseDir)
#	print(cls_vec["gm"])
# =============================================================================
#
civet.computeNativeTissueVolumes <- function(scanID, baseDir, civetVersion="1.1.9") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)

	# get vector of tissue volumes in stx space first
	cls_vec <- civet.computeStxTissueVolumes(scanID, baseDir)

	# compute the global rescaling factor
	scaling_factor <- civet.computeNativeToStxRescalingFactor(scanID, baseDir)
	
	# divide the stx volumes by the scaling factor
	cls_vec <- cls_vec / scaling_factor
	
	# return
	return(cls_vec)
}
#

# =============================================================================
# Purpose: 
# Create a brainView file to display cortical measures
#
# Example: 
#	dataFile : R Variable or .txt file
#	AtlasFile: AAL.csv
# 	leftAtlasVertices:  AAL_atlas_left.txt
#	rightAtlasVertices: AAL_atlas_right.txt
# =============================================================================
civet.CreateBrainViewFile = function(dataFile,atlasFile,atlasVertices,outputFileName,civetVersion="1.1.12") {

vertices = read.csv(atlasVertices,header=FALSE)
AALAtlas = read.csv(atlasFile)
reducedVertices = vertices[1:40962,1]
roiObj =  rep(0,length(reducedVertices))

#Input is a string specifying a file
if(is.character(dataFile))

{
	dataTable = read.table(dataFile)
	for (j in 1:dim(dataTable)[1])
	{
		reducedVerticesIndices = which(reducedVertices == dataTable[j,1],arr.ind=FALSE)
		roiObj[reducedVerticesIndices] = dataTable[j,2]
	}
}

#Input is a variablepmatch(names(dataFile[j],AALAtlas[,3])
else
{
	for (j in 1:length(dataFile))
	{
		atlasIndex = pmatch(names(dataFile[j]),AALAtlas[,3])
		reducedVerticesIndices = which(reducedVertices == AALAtlas[atlasIndex,1],arr.ind=FALSE)
		roiObj[reducedVerticesIndices] = dataFile[j]
	}
}
write.table(roiObj,outputFileName,FALSE,TRUE," ","\n","NA",".",FALSE,FALSE)
}

# =============================================================================
#
# Deprecated .. Deprecated .. Deprecated .. Deprecated .. Deprecated .. 
#
# =============================================================================
#

civetFilenames <- function(gf, idvar, prefix, basedir) {
	
	warning("This function has been deprecated. Pls use function civet.getAllFilenames() instead.", call.=TRUE)
	
  # designed for use with CIVET 1.1.9
  ids <- gf[,idvar]
  b <- paste(basedir, "/", ids, "/", sep="")
  gf$tissue <- paste(b, "classify/", prefix, "_", ids, "_cls_volumes.dat",
                     sep="")
  gf$structures <- paste(b, "segment/", prefix, "_", ids, "_masked.dat",
                         sep="")
  gf$left.thickness <- paste(b, "thickness/", prefix, "_", ids,
                             "_native_rms_rsl_tlink_20mm_left.txt", sep="")
  gf$right.thickness <- paste(b, "thickness/", prefix, "_", ids,
                             "_native_rms_rsl_tlink_20mm_right.txt", sep="")
  gf$rightROIthickness <- paste(b, "segment/", prefix, "_", ids,
                                "_lobe_thickness_tlink_20mm_right.dat", sep="")
  gf$leftROIthickness <- paste(b, "segment/", prefix, "_", ids,
                               "_lobe_thickness_tlink_20mm_left.dat", sep="")
  gf$rightROIarea <- paste(b, "segment/", prefix, "_", ids,
                           "_lobe_areas_right.dat", sep="")
  gf$leftROIarea <- paste(b, "segment/", prefix, "_", ids,
                          "_lobe_areas_left.dat", sep="")
  gf$GMVBM <- paste(b, "VBM/", prefix, "_", ids,
                    "_smooth_8mm_gm.mnc", sep="")
  gf$WMVBM <- paste(b, "VBM/", prefix, "_", ids,
                    "_smooth_8mm_wm.mnc", sep="")
  gf$CSFVBM <- paste(b, "VBM/", prefix, "_", ids,
                     "_smooth_8mm_csf.mnc", sep="")
  gf$GMVBMsym <- paste(b, "VBM/", prefix, "_", ids,
                    "_smooth_8mm_gm_sym.mnc", sep="")
  gf$WMVBMsym <- paste(b, "VBM/", prefix, "_", ids,
                    "_smooth_8mm_wm_sym.mnc", sep="")
  gf$CSFVBMsym <- paste(b, "VBM/", prefix, "_", ids,
                     "_smooth_8mm_csf_sym.mnc", sep="")
  
  
  return(gf)
}



civetAllROIs <- function(gf, defprefix) {

	warning("This function has been deprecated. Pls use function civet.AllROIs instead.", call.=TRUE)

  tissues <- anatGetAll(gf$tissue, NULL, method="text",
                        defs=paste(defprefix,"/","tissue-volumes.csv",sep=""),
                        drop=TRUE)
  colnames(tissues) <- paste(colnames(tissues), "volume")

  structures <- anatGetAll(gf$structures, NULL, method="text",
                           defs=paste(defprefix, "/","lobe-seg-volumes.csv",
                             sep=""),
                           drop=TRUE)
  colnames(structures) <- paste(colnames(structures), "volume")
  leftareas <- anatGetAll(gf$leftROIarea, NULL, method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="left")
  colnames(leftareas) <- paste(colnames(leftareas), "surface area")
  rightareas <- anatGetAll(gf$rightROIarea, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="right")
  colnames(rightareas) <- paste(colnames(rightareas), "surface area")
  leftthickness <- anatGetAll(gf$leftROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="left")
  colnames(leftthickness) <- paste(colnames(leftthickness), "mean thickness")
  rightthickness <- anatGetAll(gf$rightROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          drop=TRUE, side="right")
  colnames(rightthickness) <- paste(colnames(rightthickness), "mean thickness")

  vols <- data.frame(tissues, structures, leftareas, rightareas, leftthickness, rightthickness)
  return(vols)
}
