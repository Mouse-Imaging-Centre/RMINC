#' Minc False Discovery Rates
#' 
#' Takes the output of a mincLm type run and computes the False Discovery Rate on the results.
#' @name mincFDR
#' @aliases vertexFDR anatFDR
#' @param buffer The results of a mincLm type run.
#' @param columns A vector of column names. By default the threshold will
#' be computed for all columns; with this argument the computation can
#' be limited to a subset.
#' @param mask Either a filename or a numeric vector representing a mask
#' only values inside the mask will be used to compute the
#' threshold.
#' @param df The degrees of freedom - normally this can be determined
#' from the input object.
#' @param statType This should be either a "t","F","u","chisq" or "tlmer" depending upon the
#' type of statistic being thresholded.
#' @param method The method used to compute the false discovery
#' rate. Options are "FDR" and "pFDR".
#' @param ... extra parameters to pass to methods
#' @details This function uses the \code{qvalue} package to compute the
#'  False Discovery Rate threshold for the results of a \link{mincLm}
#'  computation. The False Discovery Rate represents the percentage of
#'  results expected to be a false positive. Two implementations can be
#'  used as specified by the method argument. "FDR" uses the
#'  implementation in \code{p.adjust}, whereas "pFDR" is a version of the
#'  postivie False Discovery Rate as found in John Storey's \code{qvalue}
#'  package. The main interface functions are 
#'  \itemize{
#'  \item{mincFDR.mincMultiDim}{ The workhorse function, used to compute q-values
#'  and thresholds for sets of minc volumes}
#'  \item{mincFDR.logLikRatio}{ Similar to above, but calculates thresholds by parametric
#'  bootstrap when possible}
#'  \item{mincFDR.mincSingleDim}{ Used when \link{mincLm}-like results are written out and read back in
#'  either to the same or another R session. In this case it loses it's \code{minMultiDim} class
#'  and must be converted back}
#'  \item{vertexFDR}{ Used with results of a \link{vertexLm}-like command. Results are converted internally
#'  to resemble a mincMultiDim and processed as normal}
#'  } 
#' @return A object of type \code{mincQvals} with the same number of columns
#'  as the input (or the subset specified by the columns argument to
#'  mincFDR). Each column now contains the qvalues for each voxel. Areas
#'  outside the mask (if a mask was specified) will be represented by a
#'  value of 1. The result also has an attribute called "thresholds"
#'  which contains the 1, 5, 10, 15, and 20 percent false discovery rate
#'  thresholds.
#' @seealso mincWriteVolume,mincLm,mincWilcoxon or mincTtest 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
#' # compute the False Discovery Rate
#' qvals <- mincFDR(vs)
#' }
#' @export
mincFDR <- function(buffer, ...) {
  UseMethod("mincFDR")
}


#' @export
vertexFDR <- function(buffer, method = "FDR", mask = NULL) {
  mincFDR.mincMultiDim(buffer, columns = NULL, mask = mask, df = NULL,
                       method = method)
}

#' @describeIn mincFDR mincSingleDim
#' @export
mincFDR.mincSingleDim <- function(buffer, df, mask = NULL, method = "qvalue", ...) {
  if (is.null(df)) {
    stop("Error: need to specify the degrees of freedom")
  }
  # if (length(df) == 1) {
  #   df <- c(1,df)
  # }
  
  dim(buffer) <- c(length(buffer), 1)
  mincFDR.mincMultiDim(buffer, columns=1, mask=mask, df=df, method=method, ...)
}

#' mincFDRMask
#'
#' Returns either the specified mask, the mask associated with the buffer, 
#' or a vector of ones to be used as a mask.
#' 
#' @param mask a mask file or vector to be passed to mincGetMask
#' if left null, the buffer is checked for a mask attribute, if no mask
#' is found, a vector of ones is used
#' @param buffer a buffer describing a minc volume  
#' @return a numeric mask vector
#' @export
mincFDRMask <- function(mask = NULL, buffer) {
  if (is.null(mask)) {
    mask <- attr(buffer, "mask")
  }
  cat("Using mask:", mask, "\n")
  
  # if mask is still null, create a vector of ones with length of buffer
  if (is.null(mask)) {
    mask <- vector(length=nrow(buffer)) + 1
  }
  else {
    mask <- mincGetMask(mask)
  }
  return(mask)
}

#' @describeIn mincFDR mincLogLikRatio
#' @export
mincFDR.mincLogLikRatio <- function(buffer, mask=NULL, ...) {
  cat("Computing FDR for mincLogLikRatio\n")
  df <- attr(buffer, "df")
  mask <- mincFDRMask(mask, buffer)
  ncols <- ncol(buffer)
  
  # compute the thresholds at several sig levels
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  
  # check whether the parametric bootstrap was estimated on the model
  haveParametricBootstrap <- "parametricBootstrap" %in% names(attributes(buffer))
  if (haveParametricBootstrap) {
    # we'll keep the estimate, corrected q|pvals and their upper and lower 95% conf limits
    # note: same as parametricBootstrap, at the moment has a limit of a single column of
    # chi-squared values (i.e. the comparison of two input models)
    ncols <- 4
  }
  
  # compute p and q values
  # pvals and qvals size corresponds to voxels inside mask, whereas output
  # is the size of the buffer. (ncols is changed if parametricBootstrap present)
  pvals <- matrix(nrow=sum(mask>0.5), ncol=ncols)
  output <- matrix(1, nrow=nrow(buffer), ncol=ncols)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=ncols)
  qvals <- pvals
  
  
  for (i in 1:ncols) {
    # compute qvals through R's p.adjust function
    currentDF <- df[[1]]
    if (i == 1 | haveParametricBootstrap == F) {
      currentDF <- df[[i]]
      pvals[, i] <- pchisq(buffer[mask>0.5, i], currentDF, lower.tail=F)
    }
    if (haveParametricBootstrap & i==1) {
      # the linear model computed as part of the parametric bootstrap
      # note: as with the parametricBootstrap code, at the moment the assumption is
      # that there were only two models being compared, and thus a single column of
      # chisq values in the input 
      
      pvals[,2:4] <- predict(attr(buffer, "parametricBootstrapModel"),
                             newdata=data.frame(chisq=pvals[,i]),
                             interval="confidence")
    }
    qvals[, i] <- p.adjust(pvals[, i], "fdr")
    
    tfunc <- function(x) { qchisq(max(x), currentDF, lower.tail=FALSE) }
    thresholds[,i] <- mincFDRThresholdVector(pvals[,i], qvals[,i], tfunc, p.thresholds)
    output[mask>0.5, i] <- qvals[, i]
  }
  
  if (haveParametricBootstrap) {
    columnNames <- c("Chisq approx.", "Corrected", "Corrected lwr CI", "Corrected upr CI")
  }
  else {
    columnNames <- colnames(buffer)
  }
  
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columnNames
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columnNames
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(output)
}

#' a utility function to compute thresholds
#'
#' @param pvals a vector of pvalues
#' @param qvals a vector of corrected qvalues (such as returend by p.adjust)
#' @param thresholdFunc a function that returns the threshold given a vector of pvalues
#' @param p.thresholds the pvalues at which to compute the threshold
#'
#' The function should be the quantile function for the distribution being tested. For example,
#' for the chi squared distribution the function would be:
#' tfunc <- function(x) { qchisq(max(x), df[[i]], lower.tail=FALSE) }
mincFDRThresholdVector <- function(pvals, qvals, thresholdFunc=NULL,
                                   p.thresholds = c(0.01, 0.05, 0.10, 0.15, 0.20)) {
  
  thresholds <- vector("numeric", length=length(p.thresholds))
  for (j in 1:length(p.thresholds)) {
    # compute thresholds; to be honest, not quite sure what the NA checking is about
    subTholdPvalues <- pvals[qvals <= p.thresholds[j]]
    subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
    
    if ( length(subTholdPvaluesNumbers) >= 1 ) {
      thresholds[j] <- thresholdFunc(subTholdPvaluesNumbers)
      #qchisq(max(subTholdPvaluesNumbers), df[[i]], lower.tail=FALSE)
    }
    else { thresholds[j] <- NA }
  }
  return(thresholds)
}


#' @describeIn mincFDR mincLmer
#' @export
mincFDR.mincLmer <- function(buffer, mask=NULL, ...) {
  cat("In mincFDR.mincLmer\n")
  
  # if no DF set, exit with message
  df <- attr(buffer, "df")
  if (is.null(df)) {
    stop("No degrees of freedom for mincLmer object. Needs to be explicitly assigned with mincLmerEstimateDF (and read the documentation of that function to learn about the dragons that be living there!).")
  }
  else {
    warning("Here be dragons! Hypothesis testing with mixed effects models is challenging, since nobody quite knows how to correctly estimate denominator degrees of freedom.")
  }
  
  # get the mask
  mask <- mincFDRMask(mask, buffer)
  # only compute stats on tlmer columns
  tlmerColumns <- grep("tlmer", attr(buffer, "stat-type"))
  ncolsToUse <- length(tlmerColumns)
  # sanity check to ensure that number of tlmer columns matches DF
  if (ncolsToUse != length(df)) {
    stop("Mismatch between DF and number of columns")
  }
  
  # compute p and q values
  # pvals and qvals size corresponds to voxels inside mask, whereas output
  # is the size of the buffer.
  pvals <- matrix(nrow=sum(mask>0.5), ncol=ncolsToUse)
  qvals <- pvals
  output <- matrix(1, nrow=nrow(buffer), ncol=ncolsToUse)
  
  # compute the thresholds at several sig levels
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=ncolsToUse)
  for (i in 1:ncolsToUse) {
    # compute qvals through R's p.adjust function
    pvals[, i] <- pt2(buffer[mask>0.5, tlmerColumns[i]], df[[i]])
    qvals[, i] <- p.adjust(pvals[, i], "fdr")
    output[mask>0.5, i] <- qvals[, i]
    for (j in 1:length(p.thresholds)) {
      # compute thresholds; to be honest, not quite sure what the NA checking is about
      subTholdPvalues <- pvals[qvals[,i] <= p.thresholds[j], i]
      subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
      
      if ( length(subTholdPvaluesNumbers) >= 1 ) {
        thresholds[j,i] <-qt(max(subTholdPvaluesNumbers)/2, df[[i]], lower.tail=FALSE)
      }
      else { thresholds[j,i] <- NA }
    }
  }
  
  columnNames <- colnames(buffer)[tlmerColumns]
  columnNamesQ <- sub("tvalue", "qvalue", columnNames)
  
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columnNames
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columnNamesQ
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(output)
}

#' @describeIn mincFDR anatLmerModel
#' @export 
mincFDR.anatLmerMod <-
  function(buffer, ...){
    mincFDR.mincLmer(buffer)
  }

#' @describeIn mincFDR mincMultiDim
#' @export
mincFDR.mincMultiDim <- function(buffer, columns=NULL, mask=NULL, df=NULL,
                                 method="FDR", statType=NULL, ...) {
  
  if(is.null(attr(buffer, "df"))) attr(buffer, "df") <- df
  originalMincAttrs <- mincAttributes(buffer)
  stattype <- originalMincAttrs$`stat-type`
  
  # Remove coefficients from buffer
  if(!is.null(stattype)){
    for (nStat in 1:length(stattype)) {
      if(stattype[nStat] == 'beta' || stattype[nStat] == 'R-squared' || stattype[nStat] == "logLik") {
        if(!exists('indicesToRemove')) {
          indicesToRemove = nStat 
        }
        else {
          indicesToRemove = c(indicesToRemove,nStat) 
        }
      }
    }
    if(exists('indicesToRemove')) {
      
      updatedAttrs <- originalMincAttrs
      updatedAttrs$`stat-type` <- updatedAttrs$`stat-type`[-indicesToRemove]
      updatedAttrs$dimnames[[2]] <- updatedAttrs$dimnames[[2]][-indicesToRemove]
      
      buffer <- buffer[,-indicesToRemove]
      buffer <- setMincAttributes(buffer, updatedAttrs)
    }
  }
  
  
  # must know the type of statistic we are dealing with
  knownStats <- c("t", "F", "u", "chisq", "tlmer")
  if (is.null(statType)) {
    # stat type not specified - must be an attribute to the buffer
    if (is.null(attr(buffer, "stat-type"))) {
      stop("Error: need to specify the type of statistic.")
    }
    else {
      statType <- attr(buffer, "stat-type")
    }
    # make sure that the stat type is recognized
    if (! all(statType %in% knownStats)) {
      stop("Error: not all the stat types are recognized. Currently allowed are: ",
           paste(knownStats, collapse=" "))
    }
    # make sure that there are either just one stat type
    # or as many as there are columns
    if (length(statType) == 1 & ncol(buffer) !=1) {
      statType <- rep(statType, ncol(buffer))
    }
    else if (length(statType) == ncol(buffer)) {
      # do nothing
    }
    else {
      stop("Error: stat type needs to be either a single entry or as many entries as there are columns in the buffer")
    }
  }
  
  
  if ( any(statType %in% "u")) {
    m <- attr(buffer, "m") 
    n <- attr(buffer, "n")
  }
  else {
    # need to know the degrees of freedom
    df <- attr(buffer, "df")
    if (is.null(df)) {
      if (any(statType %in% "tlmer")) {
        stop("Error: no degrees of freedom for mincLmer object. Needs to be explicitly assigned with mincLmerEstimateDF (and read the documentation of that function to learn about the dragons that be living there!).")
      }
      else {
        stop("Error: need to specify the degrees of freedom")
      }
    }
    if (length(df) == 1 & ncol(buffer) != 1) {
      df <- rep(list(df), ncol(buffer))
    }
    else if (length(df) == ncol(buffer)) {
      # do nothing
    }
    else {
      stop("Error: df needs to be of either length 1 or the same length as number of columns in the buffer")
    }
    #df <- vector(length=2)
    #df[1] <- ncol(attributes(buffer)$model) -1
    #df[2] <- nrow(attributes(buffer)$model) - ncol(attributes(buffer)$model)
  }
  
  if (is.null(columns)) {
    columns <- colnames(buffer)
    cat("\nComputing FDR threshold for all columns\n")
  }
  
  n.cols <- length(columns)
  n.row <-0
  if (is.matrix(buffer)) {
    n.row <- nrow(buffer)
  }
  else {
    n.row <- length(buffer)
  }
  
  if (is.null(mask)) {
    mask <- vector(length=n.row) + 1
  }
  else {
    mask <- mincGetMask(mask)
  }
  
  
  output <- matrix(1, nrow=n.row, ncol=n.cols)
  p.thresholds <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  thresholds <- matrix(nrow=length(p.thresholds), ncol=n.cols)
  
  if (any("tlmer" %in% statType)) {
    warning("Warning: computing p-values from a mincLmer call. Mixed-effects models are notoriously difficult to correctly obtain p-values from, so this is based on an approximation and might be incorrect. Read the documentation and, if in doubt, use log likelihood testing for a more correct approach.")
  }
  
  for (i in 1:n.cols) {
    cat("  Computing threshold for ", columns[i], "\n")
    pvals <- 0
    qobj <- vector("list", length(pvals))
    
    # convert statistics to p-values
    if (statType[i] %in% c("t", "tlmer")) {
      if (is.matrix(buffer)) {
        pvals <- pt2(buffer[mask>0.5, i], df[[i]])
      }
      
      else {
        pvals <- pt2(buffer[mask>0.5], df[[i]])
      }
    }
    else if (statType[i] == "F") {
      if (is.matrix(buffer)) {
        pvals <- pf(buffer[mask>0.5, i], df[[i]][1], df[[i]][2],
                    lower.tail=FALSE)
      }
      
      
      else {
        pvals <- pf(buffer[mask>0.5], df[[i]][1], df[[i]][2], lower.tail=FALSE)
      }
      
    }
    else if (statType[i] == "u") {
      pvals <- 1 - pwilcox(buffer[mask>0.5,i],m,n,lower.tail = FALSE)
    }
    else if (statType[i] == "chisq") {
      if (is.matrix(buffer)) {
        pvals <- pchisq(buffer[mask>0.5, i], df[[i]], lower.tail=F)
      }
      else {
        pvals <- pchisq(buffer[mask>0.5], df[[i]], lower.tail=F)
      }
    }  
    
    # determine corresponding q values
    if (method=="qvalue") {
      qobj <- qvalue::qvalue(pvals)
    }
    else if (method == "FDR" | method == "p.adjust") {
      qobj$pvalue <- pvals
      qobj$qvalue <- p.adjust(pvals, "fdr")
    }
    else if (method == "pFDR" | method == "fastqvalue") {
      qobj <- fast.qvalue(pvals)
    }
    # calculate thresholds at different sig levels
    for (j in 1:length(p.thresholds)) {
      if (statType[i] == "F") {
        subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
        subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
        # cat(sprintf("Number of sub-threshold F p-values: %d\n", length(subTholdPvalues)))
        if ( length(subTholdPvaluesNumbers) >= 1 ) {
          thresholds[j,i] <- qf(max(subTholdPvaluesNumbers), df[[i]][1], df[[i]][2], lower.tail=FALSE)
        } else { thresholds[j,i] <- NA }
      }
      else if (statType[i] %in% c("t", "tlmer")) {
        subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
        subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
        #cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
        if ( length(subTholdPvaluesNumbers) >= 1 ) {
          thresholds[j,i] <-qt(max(subTholdPvaluesNumbers)/2, df[[i]], lower.tail=FALSE)
        } else { thresholds[j,i] <- NA }
      }
      else if (statType[i] == "u") {
        subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
        subTholdPvaluesNumbers = subTholdPvalues[which(!is.na(subTholdPvalues))];
        #cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
        if ( length(subTholdPvaluesNumbers) >= 1 ) {
          thresholds[j,i] <-qwilcox(max(subTholdPvaluesNumbers),m,n,lower.tail = TRUE)
        } else { thresholds[j,i] <- NA }
      }
      
      else if (statType[i] == "chisq") {
        subTholdPvalues <- qobj$pvalue[qobj$qvalue <= p.thresholds[j]]
        #cat(sprintf("Number of sub-threshold t p-values: %d\n", length(subTholdPvalues)))
        if ( length(subTholdPvalues) >= 1 ) {
          thresholds[j,i] <-qchisq(max(subTholdPvalues), df[[i]], lower.tail=FALSE)
        } else { thresholds[j,i] <- NA }       
      }
    }
    output[mask>0.5,i] <- qobj$qvalue
  }
  
  rownames(thresholds) <- p.thresholds
  colnames(thresholds) <- columns
  attr(output, "thresholds") <- thresholds
  colnames(output) <- columns
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  attr(output, "DF") <- df
  class(output) <- c("mincQvals", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(output)
}

#' Get FDR Thresholds
#' 
#' @param qvals A \code{mincQvals} object, typically computed with \code{mincFDR}
#' methods
#' @return A matrix of thresholds, accessible with standard matrix indexing
#' @export 
thresholds <- 
  function(qvals){
    if(!inherits(qvals, "mincQvals")) stop("Input is not a mincQvals object")
    
    attr(qvals, "thresholds")
  }