anatGetFile <- function(filename, 
                        atlas, 
                        method = "jacobians", 
                        defs = getOption("RMINC_LABEL_DEFINITIONS"), 
                        dropLabels=FALSE, 
                        side = "both" ) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  out <- NULL
  tmpfile <- tempfile(pattern="RMINC-", fileext=".txt")
  if (method == "jacobians") {
    system(paste("label_volumes_from_jacobians", atlas, filename, ">", tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "labels") {
    # filename here should be a set of labels unique to this brain
    system(paste("volumes_from_labels_only", filename, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "means") {
    system(paste("compute_values_across_segmentation", "-m",
                 filename, atlas, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "sums") {
    system(paste("compute_counts_for_labels",
                  atlas, filename, ">", tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "slow_sums") {
    system(paste("compute_values_across_segmentation", "-s",
                 filename, atlas, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "text") {
    # values are already extracted and stored in a text file
    out <- read.table(filename, header=FALSE)
  }
  else {
    # unrecognized option...
    stop("Unrecognized option used for \"method\" (anatGetFile/anatGetAll). Available options are: jacobians, labels, means, sums, text.")
  }
  #cat("FILENAME:", filename, "\n")
  if (dropLabels == TRUE) {
    labels <- read.csv(defs)
    usedlabels <- 0
    if (side == "right") {
      usedlabels <- labels$right.label
    }
    else if (side == "left") {
      usedlabels <- labels$left.label
    }
    else if (side == "both") {
      usedlabels <- c(labels$left.label, labels$right.label)
    }
    out <- out[out$V1 %in% usedlabels,]
  }
  on.exit(unlink(tmpfile))
  return(out)
}

anatRenameRows <- function(anat, defs=getOption("RMINC_LABEL_DEFINITIONS")) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  defs <- read.csv(defs)
  rn <- rownames(anat)
  on <- as.character(rn)
  for (i in 1:length(rn)) {
    # if structure exists in both hemisphere with same number
    if (rn[i] %in% defs$right.label & rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste(defs$Structure[index])
      #cat("Match", rn[i], "to both:", on[i], "\n")
    }
    # left side only
    else if (rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$left.label)
      on[i] <- paste("left", defs$Structure[index])
      #cat("Match", rn[i], "to left:", on[i], "\n")
    }
    else if (rn[i] %in% defs$right.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste("right", defs$Structure[index])
      #cat("Match", rn[i], "to right:", on[i], "\n")
    }
    else {
      #cat("Something weird", rn[i], "\n")
    }
  }
  rownames(anat) <- on
  # make it an anat class
  class(anat) <- c("anatUnilateral", "anatMatrix", "matrix")
  # make sure we can map label numbers for later use
  attr(anat, "anatIDs") <- rn
  return(anat)
}

# both unilateral and bilateral matrices can be printed the same way
print.anatMatrix <- function(x, ...) {
  print.table(x)
}
###########################################################################################
#' @description Computes volumes, means, sums, and similar values across a
#' segmented atlas
#' @name anatGetAll
#' @title Get values given a set of files and an atlas
#' @param filenames A vector of filenames (strings) which contain the
#' information to be extracted at every structure in the atlas.
#' @param atlas A single filename containing the atlas definitions. This MINC
#' volume has to be of the same sampling (sizes and dimension
#' order) as the filenames specified in the first argument and
#' use a separate integer for each atlas label.
#' @param method A string specifying the way information is to be computed at
#' every voxel. See the details section for the possible options
#' and what they mean.
#' @param defs A string pointing to the filename containing the label
#' definitions. Used to map the integers in the atlas to a
#' proper name for the structure and contains additional
#' information for laterality of each structure. See \link{voxel_atlas_defs}
#' for details.
#' @param dropLabels Whether to return a value for every structure in the defs
#' or just for the ones actually contained in each file.
#' @param side Three choices, "right", "left", and "both" (the default)
#' which specify what labels to obtain.
#' @details 
#' anatGetAll needs a set of files along with an atlas and a set of
#' atlas definitions. In the end it will produce one value per label
#' in the atlas for each of the input files. How that value is
#' computed depends on the methods argument:
#' \itemize{
#'   \item{jacobians -}{ Each file contains log jacobians, and the volume for
#'   each atlas label is computed by multiplying the jacobian with
#'   the voxel volume at each voxel.
#'   }
#'   \item{labels -}{ Each file contains integer labels (i.e. same as the atlas).
#'   The volume is computed by counting the number of voxels with
#'   each label and multiplying by the voxel volume.
#'   }
#'   \item{means -}{ Each file contains an arbitrary number and the mean of all
#'   voxels inside each label is computed.
#'   }
#'   \item{sums -}{ Each file contains an aribtrary number and the sum of all
#'   voxels inside each label is computed.
#'   }
#'   \item{text -}{ Each file is a comma separated values text file and is simply
#'   read in.
#'   }
#' }
#' @return A matrix with ncols equal to the number of labels in the atlas and
#' nrows equal to the number of files.

#' @seealso anatLm,anatCombineStructures
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
#' volumes <- anatGetAll(filenames=filenames$absolute_jacobian, 
#'                       atlas="/tmp/rminctestdata/test_segmentation.mnc", 
#'                       method="jacobians",
#'                       defs="/tmp/rminctestdata/test_defs.csv")
#'}
#'@export
anatGetAll <- function(filenames, atlas, method="jacobians", 
                       defs = getOption("RMINC_LABEL_DEFINITIONS"), 
                       dropLabels = TRUE, 
                       side="both") {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  # Get output dimensions from full set of label definitions
  labeldefs <- read.csv(defs) 
  labels <- c(labeldefs$right.label, labeldefs$left.label)
  labels.sorted <- sort(labels)
  usedlabels <- labels.sorted[!duplicated(labels.sorted)]
  
  output <- matrix(nrow = length(usedlabels), ncol = length(filenames))
  rownames(output) <- usedlabels

  for (i in 1:length(filenames)){
    vol <- anatGetFile(filenames[i], atlas, method, defs, dropLabels, side)

    ### some quick error checking ###
    # get the list of all the labels that was found in the current file
    labels_found <- vol$V1
    # find labels that occurred in the file, but not in the mapping
    diff_between_file_and_mapping <- setdiff(labels_found, usedlabels)
    if(length(diff_between_file_and_mapping) > 0) {
      print(paste("WARNING: Labels found in the inputfile, but not in the mapping: ", diff_between_file_and_mapping))
    }
    # find labels that occurred in the mapping, but not in the file
    diff_between_mapping_and_file <- setdiff(usedlabels,labels_found)
    if(length(diff_between_mapping_and_file) > 0) {
      print(paste("WARNING: Labels found in the mapping, but not in the file: ", diff_between_mapping_and_file))
    }
    ### end of error checking ###

    for (j in 1:length(usedlabels)){
       if(usedlabels[j] == vol[j,1]){
          output[j,i] <- vol[j,2]
       }
       else {
          #If labels don't exist for a particular volume, set their value to zero
          toadd <- list(usedlabels[j],0)
          dim = nrow(vol)
          vol <- rbind(vol[1:j-1,], toadd, vol[j:dim,]) 
          output[j,i] <- 0
       }			
    }
  }
 
  attr(output, "atlas") <- atlas
  attr(output, "input") <- filenames
  if (! is.null(defs)) {
    output <- anatRenameRows(output, defs)
  }
  return(t(output))
}
###########################################################################################
#' @description Combines left and right labels from volumes obtained from anatGetAll call
#' @name anatCombineStructures
#' @title Combine left and right volumes
#' @param vols Matrix output from call to anatGetAll
#' @param method A string specifying the way information was computed at
#' every voxel ("jacobians","labels","means","sums")
#' @param defs A string pointing to the filename containing the label
#' definitions. Used to map the integers in the atlas to a
#' proper name for the structure and contains additional
#' information for laterality of each structure.
#' @details anatCombineStructures collapses left and right volume information into one
#' measure. If "jacobians","sums",or "labels" is selected then the sum of the left and right is produced, otherwise
#' the mean is produced.
#' @return A matrix with ncols equal to the number of collapsed labels

#' @seealso anatLm,anatGetAll
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
#' volumes <- anatGetAll(filenames=filenames$absolute_jacobian, 
#'                       atlas="/tmp/rminctestdata/test_segmentation.mnc", 
#'                       method="jacobians",
#'                       defs="/tmp/rminctestdata/test_defs.csv")
#' volumes_combined <- 
#'      anatCombineStructures(vols=volumes, 
#'                            method="jacobians",
#'                            defs="/tmp/rminctestdata/test_defs.csv")
#' }
#' @export
anatCombineStructures <- function(vols, method = "jacobians", 
                                  defs = getOption("RMINC_LABEL_DEFINITIONS")) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  labels <- read.csv(defs)
  combined.labels <- matrix(nrow=nrow(vols), ncol=nrow(labels))
  labelNumbers <- attr(vols, "anatIDs")
  for (i in 1:nrow(labels)) {
    if (labels$right.label[i] == labels$left.label[i]) {
      combined.labels[,i] <- vols[, labelNumbers == as.character(labels$right.label[i])]
    }
    else {
      if (method == "jacobians" || method == "labels" || method == "sums") {
        combined.labels[,i] <-
          vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]
      }
      else if (method == "means"){
        combined.labels[,i] <-
          (vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]) /2
      }
        
    }
  }
  colnames(combined.labels) <- labels$Structure
  class(combined.labels) <- c("anatCombined", "anatMatrix", "matrix")
  attr(combined.labels, "atlas") <- attr(vols, "atlas")
  attr(combined.labels, "definitions") <- defs
  return(combined.labels)
}
###########################################################################################
#' @description This function is used to compute an arbitrary function of every region in an anat structure.
#' @name anatApply
#' @title Apply function over anat structure
#' @param vols anatomy volumes
#' @param grouping grouping with which to perform operations
#' @param method The function which to apply [default mean]
#' @return  out: The output will be a single vector containing as many
#'          elements as there are regions in the input variable by the number of groupings
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- anatApply(gf$lobeThickness,gf$Primary.Diagnosis)
#' }
#' @export
anatApply <- function(vols, grouping, method=mean) {
  ngroups <- length(levels(grouping))
  output <- matrix(nrow=ncol(vols), ncol=ngroups)

  for (i in 1:ncol(vols)) {
    output[i,] <- tapply(vols[,i], grouping, method)
  }
  colnames(output) <- levels(grouping)
  rownames(output) <- colnames(vols)
  return(output)
}
###########################################################################################
#' Calculates statistics and coefficients for linear model of specified anat structure
#' @param formula a model formula. The RHS of the formula may contain one term with a matrix. If
#' so only the + operator may be used, and only two terms may appear on the RHS
#' @param data a data.frame containing variables in formula 
#' @param anat an array of atlas labels vs subject data
#' @param subset rows to be used, by default all are used
#' @return Returns an object containing the R-Squared,value,coefficients,F 
#' and t statistcs that can be passed directly into anatFDR. Additionally
#' has the attributes for model,stat type and degrees of freedom.
#' @seealso mincLm,anatLm,anatFDR 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' rmincLm= anatLm(~ Sex,gf,gf$lobeThickness)
#' } 
#' @export
anatLm <- function(formula, data, anat, subset=NULL) {
  
  #INITIALIZATION
  matrixFound = FALSE
  mmatrix =  matrix()
  matrixName = '';

  # Build model.frame
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  # Add a row to keep track of what subsetting does
  mf$rowcount <- seq(1, nrow(data))
  m <- match(c("formula", "data", "subset", "rowcount"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

   if(length(grep("\\$",formula[[2]])) > 0) {
	stop("$ Not Permitted in Formula")  
  }

  # note - the formula parsing appears to be implemented again 
  # in parseLmFormula; should reconcile.
  # Only 1 Term on the RHS
  if(length(formula[[2]]) == 1) {
    
    term <- data[[as.character(formula[[2]])]]
    
    if (is.matrix(term)) {
      # Save term name for later
      rows = c('Intercept',formula[[2]])
      matrixName = formula[[2]]
      matrixFound = TRUE
      data.matrix.left <- t(anat[mf[,"(rowcount)"],])
      #data.matrix.left <- vertexTable(filenames)
      data.matrix.right <- t(term)
    }  
  }
  # Multiple Terms on RHS
  else {
    for (nTerm in 2:length(formula[[2]])){
      rCommand = paste("term <- data$",formula[[2]][[nTerm]],sep="")
      if(!all(as.character(formula[[2]][[nTerm]]) %in% names(data))) {
        next
      }
      eval(parse(text=rCommand))	
      if (is.matrix(term)) {
        matrixName = formula[[2]][[nTerm]]
        matrixFound = TRUE
        #data.matrix.left <- vertexTable(filenames)
        data.matrix.left <- t(anat[mf[,"(rowcount)"],])
        data.matrix.right <- t(term)
        
      }  else  {
        tmpFormula = formula
        rCommand = paste("formula <-",formula[[1]],"~",formula[[2]][[nTerm]],sep="")
        eval(parse(text=rCommand))	
        mmatrix <- model.matrix(formula, mf)	
        formula = tmpFormula	
      }
    }
    rows = colnames(mmatrix)
    rows = append(rows,matrixName)
  }	
  
  
  
  # Call subroutine based on whether matrix was found
  if(!matrixFound) {    	
    mmatrix <- model.matrix(formula, mf)	
    data.matrix.left <- t(anat[mf[,"(rowcount)"],])
    rows = colnames(mmatrix) 
    data.matrix.right = matrix()
  }
  result <- .Call("vertex_lm_loop",data.matrix.left,data.matrix.right,mmatrix,PACKAGE="RMINC") 
  rownames(result) <- colnames(anat)
  
  # the order of return values is:
  #
  # f-statistic
  # r-squared
  # betas
  # t-stats
  #
  betaNames = paste('beta-',rows, sep='')
  tnames = paste('tvalue-',rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("anatModel", "matrix")
  attr(result, "atlas") <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2 + 1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

#' Anatomy Linear Mixed Effects Model
#' 
#' Apply a linear mixed effects model to the results of
#' \link{anatGetAll}.
#' 
#' @param formula the model formula
#' @param data the predictor variables
#' @param anat a subject by label matrix of anatomical 
#' summaries typically produced by \link{anatGetAll}
#' @param subset rows to subset
#' @return an \code{anatModel} object of statistical results
anatLmer <-
  function(formula, data, anat, subset = NULL, REML = TRUE, 
           control = lmerControl(), verbose = FALSE, start = NULL){
    mc <- mcout <- match.call()
    
    # Allow the user to omit the LHS by inserting a normal random variable
    RMINC_DUMMY_LHS <- rnorm(nrow(data))
    new_formula <- update(formula, RMINC_DUMMY_LHS ~ .)
    environment(new_formula) <- environment()
    mc$formula <- new_formula
    
    mc[[1]] <- quote(lme4::lFormula)
    
    # remove lme4 unknown arguments, since lmer does not know about them and keeping them
    # generates obscure warning messages
    mc <- mc[!names(mc) %in% c("anat", "subset")]
    
    lmod <- eval(mc)
    
    # code ripped from lme4:::mkLmerDevFun
    rho <- new.env(parent = parent.env(environment()))
    rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta", 
                                                    "Lambdat", "Lind")],
                                      n = nrow(lmod$X), list(X = lmod$X)))
    REMLpass <- if (REML) 
      ncol(lmod$X)
    else 0L
    
    
    mincLmerList <- list(lmod, mcout, control, start, verbose, rho, REMLpass)
    
    out <- t(apply(anat, 2, RMINC:::mincLmerOptimizeAndExtract, mincLmerList = mincLmerList))
    
    out[is.infinite(out)] <- 0            #zero out infinite values produced by vcov
    
    termnames <- colnames(lmod$X)
    betaNames <- paste("beta-", termnames, sep="")
    tnames <- paste("tvalue-", termnames, sep="")
    colnames(out) <- c(betaNames, tnames, "logLik", "converged")
    
    # generate some random numbers for a single fit in order to extract some extra info
    mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])), mincLmerList)
    
    attr(out, "stat-type") <- c(rep("beta", length(betaNames)), rep("tlmer", length(tnames)),
                                "logLik", "converged")
    # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
    attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
    attr(out, "REML") <- REML
    attr(out, "mincLmerList") <- mincLmerList
    attr(out, "atlas") <- attr(anat, "atlas")
    attr(out, "definitions") <- attr(anat, "definitions")
    attr(out, "anat") <- anat
    
    class(out) <- c("anatLmerModel", "anatModel", "matrix")
    
    return(out)
  }

#' Estimate Degrees of freedom for an anatLmer model
#' 
#' Estimate the degrees of freedom for an \code{anatLmerModel} object
#' produced by \link{anatLmer}. See \link{anatLmer} for more details
#' on degrees of freedom estimation for linear mixed effects models
#' @param model an \code{anatLmerModel}
#' @return the same model, now with degrees of freedom set
#' @export
anatLmerEstimateDF <- function(model) {
  # set the DF based on the Satterthwaite approximation
  
  mincLmerList <- attr(model, "mincLmerList")
  
  # estimated DF depends on the input data. Rather than estimate separately at every structure,
  # instead select a small number of structures and estimate DF for those structures, then keep the
  # min
  nstructures <- min(50, nrow(model))
  rstructures <- sample(1:nrow(model), nstructures)
  dfs <- matrix(nrow = nstructures, ncol = sum(attr(model, "stat-type") %in% "tlmer"))
  
  for (i in 1:nstructures) {
    structureData <- attr(model, "anat")[,rstructures[i]]
    
    ## It seems LmerTest cannot compute the deviance function for mincLmers
    ## in the current version, instead extract the model components from
    ## the mincLmerList and re-fit the lmers directly at each structure,
    ## Slower but yeilds the correct result
    lmod <- mincLmerList[[1]]
    lmod$fr[,1] <- structureData
    
    # Rebuild the environment of the formula, otherwise updating does not
    # work in the lmerTest code
    environment(lmod$formula)$lmod <- lmod
    environment(lmod$formula)$mincLmerList <- mincLmerList
    
    mmod <-
      lmerTest::lmer(lmod$formula, data = lmod$fr, REML = lmod$REML,
                     start = mincLmerList[[4]], control = mincLmerList[[3]],
                     verbose = mincLmerList[[5]])
    
    dfs[i,] <- 
      lmerTest::summary(mmod)$coefficients[,"df"]
  }
  
  df <- apply(dfs, 2, median)
  cat("Mean df: ", apply(dfs, 2, mean), "\n")
  cat("Median df: ", apply(dfs, 2, median), "\n")
  cat("Min df: ", apply(dfs, 2, min), "\n")
  cat("Max df: ", apply(dfs, 2, max), "\n")
  cat("Sd df: ", apply(dfs, 2, sd), "\n")
  
  attr(model, "df") <- df
  
  return(model)
}

#' Performs ANOVA on each region specified 
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param anat  an array of atlas labels vs subject data
#' @param subset rows to be used, by default all are used
#' @return Returns an array with the F-statistic for each model specified 
#' by formula with the following attributes
#' \itemize{
#' \item{model}{ design matrix}
#' \item{stat-type}{ type of statistic used}
#' \item{df}{ degrees of freedom of each statistic}
#' } 
#' @seealso mincAnova,vertexAnova 
#' @examples
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET",
#'                            TRUE,"1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' rmincAnova = anatAnova(~ Sex,gf,gf$lobeThickness);
#' } 
#' @export
anatAnova <- function(formula, data=NULL, anat=NULL, subset=NULL) {
  # Create Model
  m  <- match.call()
  mf <- match.call(expand.dots=FALSE)
  # in order to keep track of which subjects we need (when subsetting),
  # we will store the row numbers
  mf$rowcount <- seq(1, nrow(data))
  m  <- match(c("formula", "data", "subset", "rowcount"), names(mf), 0)  
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mmatrix <- model.matrix(formula, mf)
  
  # Get the data from the anatomy matrix using the same
  # subset as specified in the formula (using the rownumbers)
  anatmatrix <- t(anat[mf[,"(rowcount)"],])
  
  # same stats as used for the vertices
  result <- .Call("vertex_anova_loop", anatmatrix, mmatrix,attr(mmatrix, "assign"), PACKAGE="RMINC");
  
  rownames(result) <- colnames(anat)
  # unlike in the anatLm function, the column names here should be the terms of the model
  colnames(result) <- attr(terms(formula), "term.labels")
  class(result) <- c("anatModel", "matrix")
  attr(result, "atlas") <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "stat-type") <-  rep("F", ncol(result))
  
  # get the structure in order to get the degrees of freedom
  firstStructure <- anatmatrix[1,] 
  l <- lm.fit(mmatrix, firstStructure)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
    dflist[[i]] <- c(dfs[[i+1]], dfr)
  }
  attr(result, "df") <- dflist
  
  return(result)
}

#' Create a Statistic Volume
#' 
#' Convert an anatomy volume into a statistic volume
#' with the results of an \link{anatLm}
#' 
#' @param anat A minc volume object with a \code{labels} attribute
#' @param filename the filename for the new minc volume
#' @param column which column of the \link{anatLm} results to use
#' @return Invisibly returns the created volume
#' @export
anatCreateVolume <- function(anat, filename, column=1) {
  labels <- read.csv(attr(anat, "definitions"))
  volume <- mincGetVolume(attr(anat, "atlas"))
  newvolume <- volume
  for (i in 1:nrow(labels)) {
    newvolume[volume < labels$right.label[i] + 0.5 &
              volume > labels$right.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
    newvolume[volume < labels$left.label[i] + 0.5 &
              volume > labels$left.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
  }
  mincWriteVolume(newvolume, filename, attr(anat, "atlas"))
  
  return(invisible(newvolume))
}

#' @export
anatFDR <- function(buffer, method="FDR") {
  vertexFDR(buffer, method)
}

#' anatSummaries
#' 
#' These functions are used to compute the mean, standard deviation,
#'    		sum, or variance of every region in an anat structure.
#' @param anat anat structure.
#' @return  out: The output will be a single vector containing as many
#'          elements as there are regions in the input variable. 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- anatMean(gf$lobeThickness)
# }
#' @name anatSummaries
NULL

#' @describeIn anatSummaries mean
#' @export 
anatMean <- function(anat) {
   return(rowMeans(t(anat)))
}

#' @describeIn anatSummaries sum
#' @export
anatSum <- function(anat) {
   return(rowSums(t(anat)))
}

#' @describeIn anatSummaries variance
#' @export
anatVar <- function(anat) {
   return(apply(t(anat),1,var))
}

#' @describeIn anatSummaries standard deviation
#' @export
anatSd <- function(anat) {
   return(apply(t(anat),1,sd))
}
