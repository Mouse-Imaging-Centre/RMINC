#' Check for a known Civet Version
#' 
#' As different versions of Civet may reflect changes in filename or directory
#' structure, it is important that the civet.* functions be validated for each
#' new version of Civet.  This function checks whether a given Civet version
#' number has been validated for use by these routines.
#' 
#' If the passed Civet version cannot be validated, a warning message is sent
#' to std output.
#' 
#' @param civetVersion A string specifying the Civet version number, e.g.,
#' \dQuote{1.1.7}
#' @return Nothing is returned.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
#' @examples
#' 
#' \dontrun{
#' civetVersion <- "1.1.8"
#' civet.checkVersion(civetVersion)
#' }
#' @export
civet.checkVersion <- function(civetVersion) {
	if (! civetVersion %in% c("1.1.12", "1.1.9", "1.1.7")) {
		warning(sprintf("This function has not been tested with Civet version %s. Use at your own risk.", civetVersion), immediate.=TRUE)
	}
	return
}


#' Get Selected Civet Filenames
#' 
#' Returns either one or more Civet filenames, depending on file type.
#' 
#' The purpose of this function is to facilitate writing code requiring
#' manipulation of Civet products.  To this purpose, we have written a number
#' of convenience functions which, given the type of file desired and a path to
#' the Civet output directory, are able to determine and return the actual
#' filename(s).
#' 
#' @param scanID A string specifying the unique scan-id (and thus
#' sub-directory) within the Civet root output directory.
#' @param baseDir A string specifying the Civet root output directory.  This
#' directory will, in turn, contain all of the scanIDs.
#' @param civetVersion An optional string specifying the version of Civet used
#' to create the output.  This is significant since filenames and directory
#' structures may change across difference versions of Civet.
#' @param fullPath A boolean specifying whether the function is to return
#' either a fully-qualified path (TRUE) or just the filename without path
#' (FALSE).
#' @param smoothing A character code indicating the smoothing level used in
#' computing thickness, area, or volume e.g. "20mm"
#' @return Either a string or a list is returned, depending on the number of
#' filenames returned.  Specifically, a single filename is returned as a
#' string, whereas multiple filenames are returned as named lists.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
#' @examples
#' 
#' \dontrun{
#' library(RMINC)
#' 
#' # set Civet root path and scan-identifier
#' basePath <- "~/tmp/ADNI/civet/pipeOut"
#' scanID = "0221-M-AD"
#' 
#' # get the name of the aggregate tissue classification volume
#' # ... and then read it
#' classifyVolname <- civet.getFilenameClassify(scanID, basePath)
#' classifyVol <- mincIO.readVolume(classifyVolname)
#' 
#' # get the left and right gray matter surface filenames and then 
#' # ... print the names
#' gmSurfName <- civet.getFilenameGrayMatterSurfaces(scanID, basePath)
#' print(gmSurfName$left)
#' print(gmSurfName$right)
#' 
#' # get the various transformation file filenames
#' lin.xfmName <- civet.getFilenameLinearTransform(scanID, basePath)
#' print(lin.xfmName)
#' nlin.xfmNames <- civet.getFilenameNonlinearTransform(scanID, basePath)
#' print(nlin.xfmNames$xfm)		# name of the nlin xfm file
#' print(nlin.xfmNames$grid)		# name of the nlin grid file
#' }
#' @name civet.getFilename 
NULL

#' @describeIn civet.getFilename Tissue classification
#' @export
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


#' @describeIn civet.getFilename gray matter pve
#' @export
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

#' @describeIn civet.getFilename white matter pve
#' @export
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

#' @describeIn civet.getFilename csf pve
#' @export
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

#' @describeIn civet.getFilename Standard to T1 transform
#' @export
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


#' @describeIn civet.getFilename brain mask
#' @export
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

#' @describeIn civet.getFilename skull mask
#' @export
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


#' @describeIn civet.getFilename gray matter surfaces
#' @export
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

#' @describeIn civet.getFilename white matter surfaces
#' @export
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

#' @describeIn civet.getFilename mid surfaces
#' @export
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


#' @describeIn civet.getFilename cortical thickness
#' @export
civet.getFilenamesCorticalThickness <- 
  function(scanID
           , baseDir
           , civetVersion="1.1.9"
           , smoothing = "20mm"
           , fullPath=TRUE) {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	ctDir <- file.path(scanRoot, 'thickness')
	
	file_pattern <- paste0("native_rms_tlink_", smoothing, ".*\\.txt")
	
	files <- list.files(ctDir, pattern = file_pattern)
	
	# fully-qualified path requested?
	if ( fullPath ) {
		files <- file.path(ctDir, files)
	}
	
	return(list(left=files[1], right=files[2]))
}

#' @describeIn civet.getFilename cortical area
#' @export
civet.getFilenamesCorticalArea <- 
  function(scanID
           , baseDir
           , civetVersion="1.1.9"
           , smoothing = "40mm"
           , fullPath=TRUE) {
  #
  # check whether the Civet version has been tested
  civet.checkVersion(civetVersion)
  
  # get a list of matching filenames in the classify dir, and return
  baseDir <- path.expand(baseDir)
  scanRoot <- file.path(baseDir, scanID)
  ctDir <- file.path(scanRoot, 'surfaces')
  
  file_pattern <- paste0("mid_surface_rsl_(right|left)_native_area_", smoothing, "\\.txt")
  
  files <- list.files(ctDir, pattern = file_pattern)
  
  # fully-qualified path requested?
  if ( fullPath ) {
    files <- file.path(ctDir, files)
  }
  
  return(list(left=files[1], right=files[2]))
}

#' @describeIn civet.getFilename cortical volume
#' @export
civet.getFilenamesCorticalVolume <- 
  function(scanID
           , baseDir
           , civetVersion="1.1.9"
           , smoothing = "40mm"
           , fullPath=TRUE) {
  #
  # check whether the Civet version has been tested
  civet.checkVersion(civetVersion)
  
  # get a list of matching filenames in the classify dir, and return
  baseDir <- path.expand(baseDir)
  scanRoot <- file.path(baseDir, scanID)
  ctDir <- file.path(scanRoot, 'surfaces')
  
  file_pattern <- paste0("surface_rsl_(right|left)_native_volume_", smoothing, ".*\\.txt")
  
  files <- list.files(ctDir, pattern = file_pattern)
  
  # fully-qualified path requested?
  if ( fullPath ) {
    files <- file.path(ctDir, files)
  }
  
  return(list(left=files[1], right=files[2]))
}

#' @describeIn civet.getFilename surface curvature
#' @export
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

#' @describeIn civet.getFilename linear transform
#' @export
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

#' @describeIn civet.getFilename non-linear transform
#' @export
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
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' }
#' @export
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


		gf$CIVETFILES = filenames.df
    
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

#' Assemble vertex files for a CIVET run
#' 
#' Locate the vertex thickness, area, and volume files for a CIVET run
#' 
#' @inheritParams civet.getAllFilenames
#' @return A data.frame containing left and right thickness, area, and volume files.
#' @export
civet.vertexFilenames <-
  function(gf, idvar, prefix, basedir, append=TRUE, civetVersion="1.1.9"){
    ids <- getElement(gf, idvar)
    
    thickness_files <-
      lapply(ids, function(id){ 
        civet.getFilenamesCorticalThickness(id, baseDir = basedir, civetVersion = civetVersion) %>%
          as_data_frame
      }) %>% 
      bind_rows %>%
      rename_(left_thickness = "left"
              , right_thickness = "right")
    
    area_files <-
      lapply(ids, function(id){ 
        civet.getFilenamesCorticalArea(id, baseDir = basedir, civetVersion = civetVersion) %>%
          as_data_frame
      }) %>% 
      bind_rows %>%
      rename_(left_area = "left"
              , right_area = "right")
    
    volume_files <-
      lapply(ids, function(id){ 
        civet.getFilenamesCorticalVolume(id, baseDir = basedir, civetVersion = civetVersion) %>%
          as_data_frame
      }) %>% 
      bind_rows %>%
      rename_(left_volume = "left"
              , right_volume = "right")
    
    bind_cols(data_frame(ids = ids)
              , thickness_files
              , area_files
              , volume_files)
  }

#' Create a table of vertex measures
#' 
#' Read in the vertex data results of a CIVET run
#' 
#' @param vertex_files The data frame of vertex file names
#' produced by \link{civet.vertexFilenames}
#' @return A 6-element list of matrices:
#' \itemize{
#' \item{}
#' }
#' @export
civet.vertexTable <- function(vertex_files){
  columns_to_collect <- setdiff(names(vertex_files), "ids")
  
  n_vertices <- 
    getElement(vertex_files, columns_to_collect[1]) %>%
    .[!is.na(.)] %>%
    first %>%
    readLines %>%
    length
  
  read_or_NAs <- 
    function(file) 
      `if`(is.na(file)
           , rep(NA_real_, n_vertices)
           , as.numeric(readLines(file)))
    
  
  vertex_files %>%
    gather_("measure", "file", columns_to_collect) %>%
    mutate_(vertex_data = ~ lapply(file, read_or_NAs)) %>%
    arrange_(~ ids) %>%
    split(.$measure) %>%
    lapply(function(df){
      unlist(df$vertex_data) %>%
        matrix(ncol = nrow(df)) %>%
        `colnames<-`(df$ids)
    })
  
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
                        dropLabels=TRUE)
  colnames(tissues) <- paste(colnames(tissues), "volume")

  structures <- anatGetAll(gf$structures, NULL, method="text",
                           defs=paste(defprefix, "/","lobe-seg-volumes.csv",
                             sep=""),
                           dropLabels=TRUE)
  colnames(structures) <- paste(colnames(structures), "volume")
  leftareas <- anatGetAll(gf$leftROIarea, NULL, method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          dropLabels=TRUE, side="left")
  colnames(leftareas) <- paste(colnames(leftareas), "surface area")
  rightareas <- anatGetAll(gf$rightROIarea, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          dropLabels=TRUE, side="right")
  colnames(rightareas) <- paste(colnames(rightareas), "surface area")
  leftthickness <- anatGetAll(gf$leftROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          dropLabels=TRUE, side="left")
  colnames(leftthickness) <- paste(colnames(leftthickness), "mean thickness")
  rightthickness <- anatGetAll(gf$rightROIthickness, NULL,method="text",
                          defs=paste(defprefix,"/","lobe-seg-defs.csv",sep=""),
                          dropLabels=TRUE, side="right")
  colnames(rightthickness) <- paste(colnames(rightthickness), "mean thickness")

  vols <- data.frame(tissues, structures, leftareas, rightareas, leftthickness, rightthickness)
  return(vols)
}
# =============================================================================
# Purpose: 
#	Organizes CIVET .dat file based on an atlas
# Requires the atlas to be a csv file with a header and following the format
# Column 1: Numeric labels
# Column 3: Label Description
#
# Example:
#	atlasFile <- "AAL.csv"
#	dataFiles <- list of data files to organize
#
# =============================================================================
#' @title Organizes CIVET .dat files based on an Atlas
#' @description Uses an atlas to associate the measurement results in
#' the CIVET .dat output files with particular structures 
#' @param atlasFile Character path to a key to the atlas used when running civet.
#' the key should be a comma separated file with a header and the following form \cr
#' Column 1: Numeric label
#' Column 3: Corresponding structure
#' @param dataFiles character containing paths to .dat files of interest
#' typically generated with \link{civet.getAllFilenames}
#' @param civetVersion character code for the version of civet used to generate
#' the data files
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
#' @title Read all CIVET files into R
#' @usage civet.readAllCivetFiles(atlasFile, gf)
#' @param gf Data Frame containing list of all CIVET file names, and where results will be stored
#' requires gf to have an element (column) called CIVETFILES which is a data.frame containing
#' paths to the civetFiles, typically generated with \link{civet.getAllFilenames}
#' @param atlasFile Character path to a key to the atlas used when running civet.
#' the key should be a comma separated file with a header and the following form \cr
#' Column 1: Numeric label \cr
#' Column 3: Corresponding structure
#' @details Prior to running, civet.getAllFilenames may be called to generate the input argument gf .
#' This function will extract the following information from the CIVET pipeline: Lobe Area (40 mm),
#' Lobe Thickness, Lobe Volume, GI, Mid Surface Native Area, Surface Native Volume, Brain Volume Native RMS RSL tLink (20mm), Native RMS tLink (20 mm)
#' @return Returns gf augmented with additional columns containing 
#' \itemize{
#' \item{lobeArea40mm:}{ A subjects by region matrix of lobe areas parcellated by atlas region}
#' \item{lobeThickness:}{ As above, but for thicknesses}
#' \item{lobeVolume:}{ As above, but for volumes}
#' \item{GI:}{ A subjects by 6 matrix. 6 columns are left hemisphere gyrification indices for the gray matter 
#'  surface, white matter surface, midsurface of the two, followed by the same indices for the two}
#' \item{BrainVolume:}{ A subjects by 3 matrix. Three columns are CSF volume, grey matter 
#'  and white matter respectively}
#' \item{midSurfaceNativeArea:}{ A subjects by vertices matrix of mid-surface areas}
#' \item{SurfaceNativeVolume:}{ As above, but for native space volumes}
#' \item{nativeRMS_RSLtlink20mm:}{ As above, but for RMS_RSL tlink 20mm thicknesses}
#' \item{nativeRMStlink20mm:}{ As above, but for RMS tlink 20mm thicknesses}
#' }
#' @seealso  civet.getAllFilenames
#' @examples
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET", TRUE, "1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' }
#' @export
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


#' Read Civet-Generated Dat Files
#' 
#' Returns a list containing the contents of various Civet-generated text files
#' (.dat, .txt).
#' 
#' The Civet pipeline produces a number of files during its execution.  The
#' purpose of this function is to read the contents of these files and return
#' the most significant values in a list.
#' 
#' @param scanID A string specifying the unique scan-id (and thus
#' sub-directory) within the Civet root output directory.
#' @param baseDir A string specifying the Civet root output directory.  This
#' directory will, in turn, contain all of the scanIDs.
#' @param civetVersion An optional string specifying the version of Civet used
#' to create the output.  This is significant since filenames and directory
#' structures may change across difference versions of Civet.
#' @return A list is returned containing the following components:
#' \item{native_tissue_volumes}{The volumes, in cubic millimeters, of the 3
#' primary classification components. Note that the values that are returned
#' reflect the volumes within the \bold{native} brain.  These values are
#' computed by dividing the stereotactic volume values by a re-scaling factor
#' comprised of xScale * yScale * zScale (as defined in the linear xfm file).}
#' \item{gyrification_index}{The gyrification index is computed per hemisphere
#' and reflects the degree of gyrification at the cortical surface.  Value are
#' computed by dividing the cortical gray matter area by the area of a convex
#' (smooth) hull over the same area.  This computation will always yield a
#' number greater than 1, with larger numbers indicating greater gyrification.}
#' \item{native_cerebral_volumes}{While 3 values are returned, only one is a
#' particular use.  The value labeled ``cortical_gray'' reflects the volume of
#' all cortical gray matter in the native space brain.  The
#' ``extra_cerebral_csf'' value reflects the volume of all extra-cerebral csf
#' (i.e., that at the cortical surface and within the sulci), and the
#' ``wmSurface_plus_contents'' value reflects the volume of the white matter
#' surface and the volume of all components encapsulated by that surface. All
#' volume measurements are relative to the \bold{native} space brain.}
#' \item{quality_control}{Most of these values are only of interest to the
#' Civet developers.  The values labelled ``classify_qc'' reflect the
#' proportion of the various components identified by tissue classification.
#' As these values are percentages, they are applicable to both native and
#' stereotactic space.}
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
#' @examples
#' 
#' \dontrun{
#' library(RMINC)
#' 
#' # set Civet root path and scan-identifier
#' basePath <- "~/tmp/ADNI/civet/pipeOut"
#' scanID = "0221-M-AD"
#' 
#' # read in the dat files contents
#' myDats <- civet.readCivetDatFiles(scanID, basePath)
#' print(myDats)
#' 
#' # print GI info
#' print(myDats$gyrification_index)
#' print(myDats$gyrification_index["lh"])
#' 
#' # print and extract cortical volume
#' print(myDats$native_cerebral_volumes)
#' native_space_cortical_gray_volume <- myDats$native_cerebral_volumes["cortical_gray", "volume"]
#' print(native_space_cortical_gray_volume)
#' 
#' }
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

#' Compute GM, WM, and CSF Tissue Volumes
#' 
#' Returns a named vector of tissue volumes.
#' 
#' Actually, this function really returns the number of voxels of each tissue
#' type contained within the final discrete Civet-produced classification
#' volume. Now, given that Civet volumes are currently sampled using 1-mm
#' isotropic voxels, the voxel count value should also reflect the volume in
#' cubic millimeters. If this ever changes, we're going to have to make a minor
#' change in this function. Please let me know if this ever happens. The native
#' volume measurements are created by taking the stereotactic volumes and
#' dividing each of them by the xfm-derived rescaling factor.
#' 
#' @param scanID A string specifying the unique scan-id (and thus
#' sub-directory) within the Civet root output directory.
#' @param baseDir A string specifying the Civet root output directory.  This
#' directory will, in turn, contain all of the scanIDs.
#' @param civetVersion An optional string specifying the version of Civet used
#' to create the output.  This is significant since filenames and directory
#' structures may change across difference versions of Civet.
#' @return A named vector containing a value for each of the 3 tissue types.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
#' @examples
#' 
#' \dontrun{
#' library(RMINC)
#' 
#' # set Civet root path and scan-identifier
#' basePath <- "~/tmp/ADNI/civet/pipeOut"
#' scanID = "0221-M-AD"
#' 
#' # print gray matter volume in stereotactic space
#' stx_cls_vec <- civet.computeStxTissueVolumes(scanID, baseDir)
#' print(stx_cls_vec["gm"])
#' 
#' # print csf volume in native space
#' native_cls_vec <- civet.computeNativeTissueVolumes(scanID, baseDir)
#' print(native_cls_vec["csf"])
#' }
#' @name civet.computeTissueVolumes
NULL

#' @describeIn civet.computeTissueVolumes standard space
#' @export
civet.computeStxTissueVolumes <- function(scanID, baseDir, civetVersion="1.1.9") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)

	# get classify volume name
	filename <- civet.getFilenameClassify(scanID, baseDir)
	cls_vol <- mincGetVolume(filename)

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

#' @describeIn civet.computeTissueVolumes native space
#' @export
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

#' Compute Native to Stereotactic Rescaling Factor
#' 
#' Returns a single float scalar reflecting a global rescaling factor needed to
#' transform the native image to stereotactic space.
#' 
#' XFM files contain information to transform a given volume to a model. In the
#' case of Civet and rescaling, the XFM contains the rescaling factors (x,y,z)
#' needed to transform the Native volume to the model, which currently, is
#' usually the symmetrical icbm-152 model.
#' 
#' This functuon serves to compute a global rescaling factor by reading the
#' individual x,y,z rescales from the linear XFM, and returning the product.
#' 
#' Interpretation of rescaling factors: 
#' \itemize{ 
#' \item{rescale > 1.0 The native brain is \emph{expanded} to fit model} 
#' \item{rescale < 1.0 The native brain is \emph{reduced} to fit model}
#' }
#' 
#' @param scanID A string specifying the unique scan-id (and thus
#' sub-directory) within the Civet root output directory.
#' @param baseDir A string specifying the Civet root output directory.  This
#' directory will, in turn, contain all of the scanIDs.
#' @param civetVersion An optional string specifying the version of Civet used
#' to create the output.  This is significant since filenames and directory
#' structures may change across difference versions of Civet.
#' @return A scalar float reflecting the rescaling factor is returned.
#' @author Jim Nikelski \email{nikelski@@bic.mni.mcgill.ca}
#' @examples
#' \dontrun{
#' library(RMINC)
#' 
#' # set Civet root path and scan-identifier
#' basePath <- "~/tmp/ADNI/civet/pipeOut"
#' scanID = "0221-M-AD"
#' 
#' # compute the global rescaling factor
#' rescale <- civet.computeNativeToStxRescalingFactor(scanID, basePath)
#' print(rescale)
#' }
#' @export
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


#' Create a brain view file
#' 
#' Creates a text file that can be loaded into brain-view2
#' 
#' @param dataFile Either the name of a file with atlas labeling or an R array with atlas labeling
#' @param atlasFile Text file with map between atlas labels and numbers
#' @param atlasVertices Text file with map between vertex points and atlas numbers
#' @param outputFileName path to file where output will be saved 
#' @param civetVersion Version of CIVET used (Default 1.1.12)
#' @details
#' This function will create a txt file that can be loaded into brain-view2, in order to 
#' visualize the results from CIVET. This function either accepts a text file or R variable as input. 
#' If using an R variable the rows/columns must be labeled with the Atlas Names. 
#' The names are then matched to numbers as given by the AtlasFile, and then the numbers are 
#' matched to vertices in the AtlasVertexFile. 
#' @examples
#' \dontrun{
#' gf = read.csv("~/SubjectTable.csv") 
#' civet.getAllFilenames(gf,"ID","ABC123","~/CIVET", TRUE, "1.1.12") 
#' gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
#' civet.CreateBrainViewFile(gf$lobeThickness[1,],
#'                           "/Atlases/AAL/AAL.csv",
#'                           "/Atlases/AAL/AAL_atlas_left.txt",
#'                           "leftLobeThickness.txt")
#' }
#' @export
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

#' @description Create a brainview file of a specific ROI
#' @title civet.CreateBrainViewROI
#' @param atlasFile File with Atlas label (string-label(int)) mappings
#' @param atlasVertices File with vertex label(int) mappings
#' @param region Region to plot -> Must match name in atlas file exactly
#' @param civetVersion Version of CIVET 
#' @details Create a .txt file that can be overlaid unto the model template (usually in the CIVET models directory)
#' Currently only CIVET version 1.1.12 is supported.
#' @seealso civet.CreateBrainViewFile
#' @examples
#' \dontrun{
#' getRMINCTestData() 
#' civet.CreateBrainViewROI("/tmp/rminctestdata/AAL.csv",
#'                          "/tmp/rminctestdata/AAL_atlas_left.txt",
#'                          "Left Insula")
#' q()
#' #(The .txt file is written in the working directory and can be viewed via the command)
#' #brain-view2 $CIVET_DIR/models/surf_reg_model_left.obj ./LeftInsula.txt
#' }
#' @export
civet.CreateBrainViewROI <- function(atlasFile,atlasVertices,region,civetVersion="1.1.12") {
	
	if ( civetVersion != "1.1.12" ) {
		warning(sprintf("This function has not been tested with Civet version %s. Use at your own risk.", civetVersion), immediate.=TRUE)
	}

	vertices = read.csv(atlasVertices,header=FALSE)
	AALAtlas = read.csv(atlasFile)
	reducedVertices = vertices[1:40962,1]
	roiObj =  rep(0,length(reducedVertices))
	
	# Find index of ROI
	roiIndex <- which(AALAtlas[,3] == region)
	
	# Match against label
	roiLabel <- which(reducedVertices == AALAtlas[roiIndex,1],arr.ind=FALSE)
	
	# Fill in label
	roiObj[roiLabel] = 1

	# Write out (remove spaces)
	write.table(roiObj,gsub(" ","",paste(region,".txt",sep="")),FALSE,TRUE,"","\n","NA",".",FALSE,FALSE)
}

#'@title Flatten CIVET results for dplyr
#'@description
#'Convert the data.frame/Matrix/list fusion object produced by civet.readAllCivetFiles
#'to a data.frame usable with dplyr etc.
#'@param civetResults A data.frame produced by \link{civet.readAllCivetFiles}
#'@param columnsToKeep vector of column names or indices of columns from the original
#'frame to copy to the normalized table. Columns not produced by \link{civet.readAllCivetFiles}
#'and not specified here are dropped. See details
#'@details The columnsToKeep vector needs to include the subject identifier or it will be
#'dropped. Ideally other columns in the vector should be proper vector/list columns compatible
#'with dplyr
#'@return data.frame containing results of \link{civet.readAllCivetFiles}, all sub-data.frames
#'and matrices are expanded, non-standard characters in column names are replaced with underscores,
#'and a prefix denoting the origin sub-data is given.  
#'@export
civet.flattenForDplyr <-
  function(civetResults, columnsToKeep){
    
    cleanAndPrefixColNames <- 
      function(mat){
        matrixName <- deparse(substitute(mat))
        colnames(mat) %>% 
          gsub("[ [:punct:]]+", "_", .) %>%
          gsub("^_|_$", "", .) %>%
          paste(matrixName, ., sep = "_")
      }
    
    normalized_frame <- civetResults[ , columnsToKeep, drop = FALSE]
    
    lobeArea40mm <- civetResults$lobeArea40mm
    colnames(lobeArea40mm) <- cleanAndPrefixColNames(lobeArea40mm)
    lobeArea40mm <- lobeArea40mm %>% as.data.frame

    
    lobeThickness <- civetResults$lobeThickness
    colnames(lobeThickness) <- cleanAndPrefixColNames(lobeThickness)
    lobeThickness <- lobeThickness %>% as.data.frame
    
    lobeVolume <- civetResults$lobeVolume
    colnames(lobeVolume) <- cleanAndPrefixColNames(lobeVolume)
    lobeVolume <- lobeVolume %>% as.data.frame
    
    GI <- civetResults$GI
    colnames(GI) <- cleanAndPrefixColNames(GI)
    GI <- GI %>% as.data.frame
    
    BrainVolume <- civetResults$BrainVolume
    colnames(BrainVolume) <- cleanAndPrefixColNames(BrainVolume)
    BrainVolume <- BrainVolume %>% as.data.frame
    
    midSurfaceNativeArea <- civetResults$midSurfaceNativeArea
    midSurfaceNativeArea <- 
      lapply(1:nrow(midSurfaceNativeArea), function(i) midSurfaceNativeArea[i,])
    
    nativeRMS_RSLtlink20mm <- civetResults$nativeRMS_RSLtlink20mm
    nativeRMS_RSLtlink20mm <-
      lapply(1:nrow(nativeRMS_RSLtlink20mm), function(i) nativeRMS_RSLtlink20mm[i,])
    
    nativeRMStlink20mm <- civetResults$nativeRMStlink20mm
    nativeRMStlink20mm <-
      lapply(1:nrow(nativeRMStlink20mm), function(i) nativeRMStlink20mm[i,])
    
    normalized_frame <-
      normalized_frame %>%
      bind_cols(lobeArea40mm, lobeVolume, lobeThickness, GI, BrainVolume) %>%
      mutate(midSurfaceNativeArea, nativeRMS_RSLtlink20mm, nativeRMStlink20mm)
  }

#' Read CIVET QC data
#' 
#' After running the CIVET quality control pipeline, import the results
#' 
#' @param basedir The CIVET output directory
#' @param civetVersion the version of CIVET used, currently only supports
#' 1.1.12
#' @return A table of QC results including whether or not the subjects passed
#' overall quality control. See \url{http://www.bic.mni.mcgill.ca/ServicesSoftware/QualityControlCIVET12}
#' for more details.
#' @export
civet.readQC <-
  function(basedir, civetVersion="1.1.12"){
    if(civetVersion != "1.1.12")
      warning("Unsure how to deal with QC results for civet version: ", civetVersion,
              "\n trying 1.1.12 approach")
    
    qc_file <- 
      file.path(basedir, "QC") %>%
      list.files(full.names = TRUE) %>%
      grep("\\.glm", ., value = TRUE)
    
    if(length(qc_file) == 0) stop("No QC data found")
    if(length(qc_file) > 1) stop("More than one results table found, aborting")
    
    good_med_bad <- 
      function(vec, mb, bb)  
        ifelse(vec < mb, "good", ifelse(vec < bb, "med", "bad"))
    
    qc_res <- 
      read.table(qc_file, stringsAsFactors = FALSE)[,1:24] %>%
      setNames(c("ID", "x-step", "y-step", "z-step", "StxMaskErr", "CSFcls", 
                 "GMcls", "WMcls", "GMCortical", "WMsurf", "GMsurf", "SRLeft", 
                 "SRRight", "SSLeft", "SSRight", "MeanCTLeft", "MeanCTRight", 
                 "MeanWM-T1", "StdDevWM-T1", "MeanGM-T1", "StdDevGM-T1", "GIleft", 
                 "GIright", "GIfull")) %>%
      mutate_(
        CSFcls_score   = ~ ifelse(between(CSFcls, 15, 80), "good", "med"),
        GMcls_score    = ~ ifelse(between(GMcls, 15, 80), "good", "med"),
        WMcls_score    = ~ ifelse(between(WMcls, 15, 80), "good", "med"),
        WMsurf_score   = ~ good_med_bad(WMsurf, 10, 20),
        GMsurf_score   = ~ good_med_bad(GMsurf, 10, 20),
        SRLeft_score   = ~ good_med_bad(SRLeft, 250, 500),
        SRRight_score  = ~ good_med_bad(SRRight, 250, 500),
        SSLeft_score   = ~ good_med_bad(SSLeft, 250, 500),
        SSRight_score  = ~ good_med_bad(SSRight, 250, 500)
      ) %>%
      cbind(
        rowwise(.) %>% 
          do(QC_PASS = 
               as_data_frame(.) %>%
               select(matches("_score")) %>% 
               unlist %>%
               `!=`("bad") %>%
               all) %>%
          mutate_(QC_PASS = ~ unlist(QC_PASS))
      )
  }
