civetFilenames <- function(gf, idvar, prefix, basedir) {
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
