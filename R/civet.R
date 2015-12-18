
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
civet.getFilenameCorticalThickness <- function(scanID, baseDir, civetVersion="1.1.9", 
                                               fullPath=TRUE, blur = "20mm") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	ctDir <- file.path(scanRoot, 'thickness')
	files <- 
	  list.files(ctDir, 
	             pattern = glob2rx(
	               sprintf("*_native_rms_rsl_tlink_%s_*.txt", blur)))
	
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
civet.getFilenameMeanSurfaceCurvature <- function(scanID, baseDir, civetVersion="1.1.9", 
                                                  fullPath=TRUE, blur = "20mm") {
	#
	# check whether the Civet version has been tested
	civet.checkVersion(civetVersion)
	
	# get a list of matching filenames in the classify dir, and return
	baseDir <- path.expand(baseDir)
	scanRoot <- file.path(baseDir, scanID)
	ctDir <- file.path(scanRoot, 'thickness')
	files <- 
	  list.files(ctDir, 
	             pattern=glob2rx(
	               sprintf("*_native_mc_rsl_%s_*.txt", blur)))
	
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

civet.getAllFilenames <- function(gf, idvar, prefix, basedir, 
                                  append=TRUE, civetVersion="1.1.9",
                                  blurs = civet.blurConfigure(version = civetVersion)) {
  # designed for use with CIVET 1.1.9 and CIVET 1.1.12
  if ( civetVersion != "1.1.9"  && civetVersion != "1.1.12" ) {
    warning("This function has only been tested with Civet version 1.1.9. and 1.1.12 Use at your own risk.", immediate.=TRUE)
  }
  
  if (civetVersion == "1.1.12")
  {	
    return(civet.getAllFilenames_1_1_12(gf = gf, 
                                        idvar = idvar, 
                                        prefix = prefix, 
                                        basedir = basedir, 
                                        blurs = blurs))
  } else {
    return(civet.getAllFilenames_other(gf = gf,
                                       idvar = idvar,
                                       prefix = prefix,
                                       append = append,
                                       basedir = basedir,
                                       blurs = blurs))
  }
}

civet.blurConfigure <- function(..., version = "1.1.12"){
  
  #Default blur values for CIVET version 1.1.12
  default_blurs_1_1_12 <-
    list(
      lobeArea = "40mm",
      lobeThickness = "20mm",
      lobeVolume = "40mm",
      midSurfaceNativeArea = "40mm",
      surfaceNativeVolume = "40mm",
      nativeRMS_RSLtlink = "20mm",
      nativeRMStlink = "20mm"
    )
  
  #Default blur values for anything else (really just version 1.1.9)
  default_blurs_other <-
    list(
      thickness = "20mm",
      ROIthickness = "20mm",
      GMVBM = "8mm",
      WMVBM = "8mm",
      CSFVBM = "8mm"
    )
  
  if(version == "1.1.12"){
    blurs <- default_blurs_1_1_12
  } else {
    blurs <- default_blurs_other
  }
  
  #To allow users to optionally set blur variables convert the list
  #to an environment so the users ... can be assigned with list2env later
  blurs <- list2env(blurs)
  
  userOptions <- list(...)
  
  #If no ... arguments return all blurs 
  if(length(userOptions) == 0) return(as.list(blurs))

  #If ... arguments are unnamed return the specified blurs
  if(is.null(names(userOptions))) {
    return(as.list(blurs)[unlist(userOptions)])
  }
  
  #Perform name partial matching, replacing the partial matches
  #with their full names, die if unknown parameter is encountered
  #I did not use match.arg for better error messaging
  names(userOptions) <-
    vapply(names(userOptions), function(name){
      
      matched_name_index <- pmatch(name, ls(blurs))
      matched_name <- ls(blurs)[matched_name_index]
      
      if(is.na(matched_name)){
        stop(
          sprintf("%s is not in the list of configurable blurs for CIVET version %s",
                  name, version), call. = FALSE)
      }
      
      return(matched_name)
    }, character(1))
  
  #Assign user supplied parameters to the blurs environment
  list2env(userOptions, envir = blurs)
  
  #Return the new parameters as a list
  return(as.list(blurs))
}

civet.getAllFilenames_1_1_12 <- 
  function(gf, idvar, prefix, basedir, 
           blurs = civet.blurConfigure(version = "1.1.12")){
    
    # extract the scanIDs from the glim
    ids <- gf[,idvar]
    
    # create the scan-level root dir names
    scan_dir <- file.path(basedir, ids)
    surface_prefixed <- file.path(scan_dir, "surfaces", paste0(prefix, "_"))
    classify_prefixed <- file.path(scan_dir, "classify", paste0(prefix, "_"))
    thickness_prefixed <- file.path(scan_dir, "thickness", paste0(prefix, "_"))
    
    build_filename <- function(dir_prefixed, file) paste0(dir_prefixed, ids, file)

    filenames <- list()
    
    filenames <-
      within(
        filenames,
        {
          leftGIFiles  <-
            build_filename(surface_prefixed, "_gi_left.dat")
          
          rightGIFiles <-
            build_filename(surface_prefixed, "_gi_right.dat")
          
          leftlobeArea40mmFiles  <- 
            build_filename(surface_prefixed, 
                           sprintf("_lobe_areas_%s_left.dat",
                                   blurs$lobeArea))
          
          rightlobeArea40mmFiles  <- 
            build_filename(surface_prefixed,
                           sprintf("_lobe_areas_%s_right.dat",
                                   blurs$lobeArea))
          
          leftlobeThicknessFiles  <- 
            build_filename(surface_prefixed,
                           sprintf("_lobe_thickness_tlink_%s_left.dat",
                                   blurs$lobeThickness))
          
          rightlobeThicknessFiles  <-
            build_filename(surface_prefixed,
                           sprintf("_lobe_thickness_tlink_%s_right.dat",
                                   blurs$lobeThickness))
          
          leftlobeVolumeFiles  <- 
            build_filename(surface_prefixed,
                           sprintf("_lobe_volumes_%s_left.dat",
                                   blurs$lobeVolume))
          
          rightlobeVolumeFiles  <- 
            build_filename(surface_prefixed, 
                           sprintf("_lobe_volumes_%s_right.dat",
                                   blurs$lobeVolume))
          
          midSurfaceleftNativeArea <- 
            build_filename(surface_prefixed,
                           sprintf("_mid_surface_rsl_left_native_area_%s.txt",
                                   blurs$midSurfaceNativeArea))
          
          midSurfacerightNativeArea <- 
            build_filename(surface_prefixed, 
                           sprintf("_mid_surface_rsl_right_native_area_%s.txt",
                                   blurs$midSurfaceNativeArea))
          
          SurfaceleftNativeVolume <- 
            build_filename(surface_prefixed, 
                           sprintf("_surface_rsl_left_native_volume_%s.txt",
                                   blurs$surfaceNativeVolume))
          
          SurfacerightNativeVolume <- 
            build_filename(surface_prefixed, 
                           sprintf("_surface_rsl_right_native_volume_%s.txt",
                                   blurs$surfaceNativeVolume)) 
          
          brain_volume <- 
            build_filename(classify_prefixed, "_cls_volumes.dat")
          
          cerebral_volume <- 
            build_filename(thickness_prefixed, "_cerebral_volume.dat")
          
          assign(
            sprintf("nativeRMS_RSLtlink%sleft", blurs$nativeRMS_RSLtlink),
            build_filename(thickness_prefixed, 
                           sprintf("_native_rms_rsl_tlink_%s_left.txt",
                                   blurs$nativeRMS_RSL))) 
          
          assign(
            sprintf("nativeRMS_RSLtlink%sright", blurs$nativeRMS_RSLtlink),
            build_filename(thickness_prefixed, 
                           sprintf("_native_rms_rsl_tlink_%s_right.txt",
                                   blurs$nativeRMS_RSL)))  
          
          assign(
            sprintf("nativeRMStlink%sleft", blurs$nativeRMStlink), 
            build_filename(thickness_prefixed, 
                           sprintf("_native_rms_tlink_%s_left.txt",
                                   blurs$nativeRMStlink))) 
          
          assign(
            sprintf("nativeRMStlink%sright", blurs$nativeRMStlink), 
            build_filename(thickness_prefixed, 
                           sprintf("_native_rms_tlink_%s_right.txt",
                                   blurs$nativeRMStlink)))
        })
    
    #Reverse column ordering to match old getAllFilenames
    #cast as data.frame
    filenames.df <- as.data.frame(rev(filenames), stringsAsFactors = FALSE)
    
    gf$CIVETFILES <- filenames.df
    
    return(gf)
  }

civet.getAllFilenames_other <-
  function(gf, idvar, prefix, basedir, append,
           blurs = civet.blurConfigure(version = "other")){
    
    ids <- gf[,idvar]
    
    # create the scan-level root dir names
    scan_dir <- file.path(basedir, ids)
    segment_prefixed <- file.path(scan_dir, "segment", paste0(prefix, "_"))
    classify_prefixed <- file.path(scan_dir, "classify", paste0(prefix, "_"))
    thickness_prefixed <- file.path(scan_dir, "thickness", paste0(prefix, "_"))
    vbm_prefixed <- file.path(scan_dir, "VBM", paste0(prefix, "_"))
    
    build_filename <- function(dir_prefixed, file) paste0(dir_prefixed, ids, file)
    
    filenames <- list()
    
    filenames <-
      within(
        filenames,
        {
          tissue <- build_filename(classify_prefixed, "_cls_volumes.dat")
          structures <- build_filename(segment_prefixed, "_masked.dat")
          left.thickness <- 
            build_filename(thickness_prefixed,
                           sprintf("_native_rms_rsl_tlink_%s_left.txt",
                                   blurs$thickness))
          right.thickness <- 
            build_filename(thickness_prefixed,
                           sprintf("_native_rms_rsl_tlink_%s_right.txt",
                                   blurs$thickness))
          rightROIthickness <- 
            build_filename(segment_prefixed,
                           sprintf("_lobe_thickness_tlink_%s_right.dat",
                                   blurs$ROIthickness))
          leftROIthickness <- 
            build_filename(segment_prefixed,
                           sprintf("_lobe_thickness_tlink_%s_left.dat",
                                   blurs$ROIthickness))
          rightROIarea <- 
            build_filename(segment_prefixed,"_lobe_areas_right.dat")
          
          leftROIarea <- 
            build_filename(segment_prefixed, "_lobe_areas_left.dat")
          
          GMVBM <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_gm.mnc",
                                   blurs$GMVBM))
          WMVBM <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_wm.mnc",
                                   blurs$WMVBM))
          CSFVBM <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_csf.mnc",
                                   blurs$CSFVBM))
          GMVBMsym <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_gm_sym.mnc",
                                   blurs$GMVBM))
          WMVBMsym <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_wm_sym.mnc",
                                   blurs$WMVBM))
          CSFVBMsym <- 
            build_filename(vbm_prefixed,
                           sprintf("_smooth_%s_csf_sym.mnc",
                                   blurs$CSFVBM))
        }
      )
    
    filenames.df <- as.data.frame(rev(filenames), stringsAsFactors = FALSE)
    
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
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET", TRUE, "1.1.12")
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

###########################################################################################
#' @description Create a brainview file of a specific ROI
#' @name civet.CreateBrainViewROI
#' @title civet.CreateBrainViewROI
#' @param atlasFile File with Atlas label (string-label(int)) mappings
#' @param atlasVertices File with vertex label(int) mappings
#' @param region Region to plot -> Must match name in atlas file exactly
#' @param civetVersion Version of CIVET 
#' @details Create a .txt file that can be overlaid unto the model template (usually in the CIVET models directory)
#' Currently only CIVET version 1.1.12 is supported.
#' @seealso civet.CreateBrainViewFile
#' @examples
#' getRMINCTestData() 
#' civet.CreateBrainViewROI("/tmp/rminctestdata/AAL.csv","/tmp/rminctestdata/AAL_atlas_left.txt","Left Insula")
#' q()
#' (The .txt file is written in the working directory and can be viewed via the command)
#' brain-view2 $CIVET_DIR/models/surf_reg_model_left.obj ./LeftInsula.txt
###########################################################################################
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
civet.flattenForDplyr <-
  function(civetResults, columnsToKeep){
    if(!require(dplyr)) stop("Please install dplyr to use this command")
    
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
      bind_cols(lobeArea40mm, lobeThickness, GI, BrainVolume) %>%
      mutate(midSurfaceNativeArea, nativeRMS_RSLtlink20mm, nativeRMStlink20mm)
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
