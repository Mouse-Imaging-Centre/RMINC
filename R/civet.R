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

  
  return(gf)
}
