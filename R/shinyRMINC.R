launch_shinyRMINC <- function(statsoutput, anatVol, volumes=NULL, keepBetas=FALSE, plotcolumns=NULL, modelfunc=NULL) {
  # output goes into temporary directory
  tmpDir <- tempdir()
  # copy two files into temporary directory
  file.copy(from = "/home/jlerch/git/RMINC-rstudio/inst/shinyRMINC/server.R", to = tmpDir, overwrite=TRUE)
  file.copy(from = "/home/jlerch/git/RMINC-rstudio/inst/shinyRMINC/ui.R", to = tmpDir, overwrite=TRUE)
  # need to turn the statsoutput into a list
  gfs <- data.frame(filenames = attributes(statsoutput)$filenames)
  if (!is.null(plotcolumns)) {
    gfs <- cbind(gfs, plotcolumns)
  }
  statsList <- list()
  sType <- attributes(statsoutput)$`stat-type`
  cNames <- colnames(statsoutput)
  cat("Converting to mincArray - might take a second or two\n")
  for (i in 1:length(cNames)) {
    if (keepBetas || sType[i] != "beta") {
      # make symmetric for betas and tvalues
      if (sType[i] %in% c("beta", "t")) {
        symmetric <- TRUE
      }
      else {
        symmetric <- FALSE
      }
      statsList[[cNames[i]]] <- list(data=mincArray(statsoutput, cNames[i]),
                                     symmetric=symmetric,
                                     legendTitle=paste(sType[i], "value"))
    }
  }
  d <- dim(anatVol)
  m <- attributes(statsoutput)$model
  if (is.null(modelfunc)) {
    if (any(sType == "t")) {
      modelfunc <- function(x) { summary(lm(x ~ m -1)) }
    }
    else {
      modelfunc <- function(x) { anoval(lm(x ~ m -1)) }
    }
  }
  # save file into tempdir
  cat("Launching shiny\n")
  shiny:::runApp(tmpDir)
}
