#' launch a shiny based inspector
#'
#' This uses shiny to inspect the output of mincLm or mincAnova. Produces
#' various views, the ability to plot individual voxels, and get the output of
#' the stats model at that voxel.
#'
#' @param statsoutput the output of mincLm or mincAnova
#' @param anatVol the anatomical volume on which to display the statistics. Can
#'   be passed a filename, a MINC volume (from mincGetVol), or a mincArray.
#' @param volumes a matrix or data frame of volumes for plotting
#' @param keepBetas whether to include the beta coefficients
#' @param plotcolumns extra data to be used for plotting
#' @param modelfunc optional modelling function
#' @param anatLow The lower threshold value for displaying the underlying anatomy
#' @param anatHigh The upper threshold value for displaying the underlying anatomy
#' @examples
#' \dontrun{
#' vs <- mincLm(reljacobians02 ~ sex*treatment, subset(gfs, treatment != "None"))
#' anatVol <- mincArray(mincGetVolume("anatomyfile.mnc"))
#' launch_shinyRMINC(vs, anatVol, volumes=gfs$vols, 
#'                   plotcolumns=gfs[,c("sex", "Neonatal")], keepBetas=F)
#' }
#' @export
launch_shinyRMINC <- function(statsoutput, anatVol, volumes=NULL, keepBetas=FALSE, plotcolumns=NULL, modelfunc=NULL
                              , anatLow = 700, anatHigh = 1400) {
  gfs <- data.frame(filenames = attributes(statsoutput)$filenames)
  if (!is.null(plotcolumns)) {
    gfs <- cbind(gfs, plotcolumns)
  }

  # if anatVol is a filename, read it in now.
  if (is.character(anatVol)) {
    anatVol <- mincArray(mincGetVolume(anatVol))
  }

  # hokey check for whether anatVol is a mincArray or not
  if (is.null(dim(anatVol))) {
    anatVol <- mincArray(anatVol)
  }
  
  if (!is.null(volumes)) {
    gfs$vols <- volumes
  }
  
  statsList <- list()
  sType <- attributes(statsoutput)$`stat-type`
  cNames <- colnames(statsoutput)
  cat("Converting to mincArray - might take a second or two\n")
  for (i in 1:length(cNames)) {
    if (keepBetas || sType[i] != "beta") {
      # make symmetric for betas and tvalues
      if (sType[i] %in% c("beta", "t", "tlmer")) {
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
      modelfunc <- function(x) { anova(lm(x ~ m -1)) }
    }
  }
  cat("Launching shiny\n")
  shiny::runApp(system.file("shinyRMINC/", package="RMINC"))
}
