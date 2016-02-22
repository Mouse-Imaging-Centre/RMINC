launch_shinyRMINC <- function(statsoutput, anatVol, volumes=NULL, keepBetas=FALSE) {
  # output goes into temporary directory
  tmpDir <- tempdir()
  # copy two files into temporary directory
  file.copy(from = "/home/jlerch/git/RMINC-rstudio/inst/shinyRMINC/server.R", to = tmpDir, overwrite=TRUE)
  file.copy(from = "/home/jlerch/git/RMINC-rstudio/inst/shinyRMINC/ui.R", to = tmpDir, overwrite=TRUE)
  # need to turn the statsoutput into a list
  gfs <- data.frame(filenames = attributes(statsoutput)$filenames)
  statsList <- list()
  sType <- attributes(statsoutput)$`stat-type`
  cNames <- colnames(statsoutput)
  cat("Converting to mincArray - might take a second or two\n")
  for (i in 1:length(cNames)) {
    if (keepBetas || sType != "beta") {
      # make symmetric for betas and tvalues
      if (sType == "beta" || sType == "t") {
        symmetric <- TRUE
      }
      else {
        symmetric <- TRUE
      }
      statsList[[cNames[i]]] <- list(data=mincArray(statsoutput, cNames[i]),symmetric=symmetric)
    }
  }
  d <- dim(anatVol)
  # save file into tempdir
  cat("Launching shiny\n")
  shiny:::runApp(tmpDir)
}
