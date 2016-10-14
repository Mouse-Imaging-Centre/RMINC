#' Minc Voxel Summary Functions
#' 
#' Compute the mean, standard deviation, sum, or variance at every voxel across a
#' a set of MINC volumes.
#' An optional grouping variable will split the computation by group
#' rather than performing it across all volumes as is the default.
#' @param filenames Filenames of the MINC volumes across which to create the
#' descriptive statistic.
#' @param grouping Optional grouping - contains same number of elements as
#' filenames; the results will then have the descriptive
#' statistic computed separately for each group, or in the case of method = "correlation"
#' this is the variable to correlate against.
#' @param mask A mask specifying which voxels are to be included in the
#' summary.
#' @param maskval the value in the mask used to select unmasked voxels, 
#' defaults to any positive intensity from 1-99999999 internally expanded to
#' .5 - 99999999.5. If a number is specified voxels with intensities 
#' within 0.5 of the chosen value are considered selected.
#' @param method the type of summarys statistic to calculate for each voxel 
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files. If a
#'          grouping factor was specified then the output will be a
#'          matrix consisiting of as many rows as there were voxels in
#'          the files, and as many columns as there were groups.
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mm <- mincMean(gf$jacobians_0.2) 
#' ms <- mincSd(gf$jacobians_0.2)
#' mv <- mincVar(gf$jacobians_0.2,gf$Strain) 
#' ms2 <- mincSum(gf$jacobians_0.2,gf$Strain)
#' }
#' @export
mincSummary <- function(filenames, grouping=NULL, mask=NULL, method="mean", maskval=NULL) {
  mincFileCheck(filenames)
  enoughAvailableFileDescriptors(length(filenames))
  
  if (is.null(grouping)) {
    grouping <- rep(1, length(filenames))
  }
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  result <- .Call("minc2_model",
                  as.character(filenames),
                  matrix(),
                  as.double(grouping)-1,
                  NULL,
                  as.double(! is.null(mask)),
                  as.character(mask),
                  as.double(minmask),
                  as.double(maxmask),
                  NULL, NULL,
                  as.character(method), PACKAGE="RMINC")
  attr(result, "likeVolume") <- as.character(filenames[1])
  attr(result, "filenames") <- as.character(filenames)
  
  
  if (is.null(grouping)) {
    class(result) <- c("mincSingleDim", "numeric")
  }
  else {
    class(result) <- c("mincMultiDim", "matrix")
    if(!grepl("t-test",method) && !grepl("correlation",method) && !grepl("wilcoxon",method)) {
      colnames(result) <- levels(grouping)
    }
  }
  return(result)
}

#' @describeIn mincSummary mean
#' @export
mincMean <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="mean", maskval=maskval)
  return(result)
}

#' @describeIn mincSummary Variance
#' @export
mincVar <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  return(result)
}

#' @describeIn mincSummary Sum
#' @export
mincSum <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="sum", maskval=maskval)
  return(result)
}

#' @describeIn mincSummary Standard Deviation
#' @export
mincSd <- function(filenames, grouping=NULL, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="var", maskval=maskval)
  result <- sqrt(result)
  return(result)
}


#' @describeIn mincSummary Correlation
#' @export 
mincCorrelation <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="correlation", maskval=maskval)
  return(result)
}

#' Perform Arbitrary calculations on a collection of mincVolumes
#' 
#' An RCPP variant of \link{mincApply}, the primary advantage being
#' that functions of an arbitrary number of arguments can be passed
#' to mincApplyRCPP.
#' @param filenames The name of the files to apply over
#' @param fun the function to apply
#' @param ... additional parameters to fun
#' @param mask a numeric mask vector
#' @param maskval An integer specifying the value inside the mask where to
#' apply the function. If left blank (the default) then anything
#' above 0.5 will be considered inside the mask. This argument
#' only works for mincApply, not pMincApply.
#' @param filter_masked Whether or not to remove the masked values
#' from the resultant object
#' @param slab_sizes a three element numeric vector indicating the size
#' in voxels of the hyperslab to read for each file. Useful for managing
#' memory use - larger slabs are faster but require more memory. Sizes
#' must be an even factor of their respective volume dimensions.
#' @param return_indices Whether to return the voxel positions of the results
#' generally for internal use only.
#' @param collate A function to (potentially) collapse the result list
#' examples include link{unlist} and \link{simplify2array}, defaulting
#' to \link{simplify2minc} which creates an object of type \code{mincMultiDim},
#' \code{mincSingleDim}, or \code{mincList} depending on the result structure.
#' @return a list of results subject the the collate function
#' @export
mincApplyRCPP <- 
  function(filenames, 
           fun, 
           ..., 
           mask = NULL, 
           maskval = NULL,
           filter_masked = FALSE,
           slab_sizes = c(1,1,1),
           return_indices = FALSE,
           collate = simplify2minc){
    
    stopifnot(!is.null(filenames), !is.null(fun))
    enoughAvailableFileDescriptors(length(filenames))
    
    apply_fun <- 
      function(x, extra_arguments)
        do.call(match.fun(fun), c(list(x), extra_arguments))
    
    args <- list(...)
    filenames <- as.character(filenames)
    
    if(is.null(maskval)){
      minmask <- 1
      maxmask <- 99999999
    } else {
      minmask <- maxmask <- maskval
    }
    
    if(is.null(mask)){
      use_mask <- FALSE
      mask = ""
    } else {
      use_mask <- TRUE
    }
    
    masked_value <- getOption("RMINC_MASKED_VALUE")
    
    results_indexed <-
      .Call("RMINC_rcpp_minc_apply",
            filenames,
            use_mask = use_mask,
            mask = mask,
            mask_lower_val = minmask,
            mask_upper_val = maxmask,
            value_for_mask = masked_value,
            filter_masked = filter_masked,
            slab_sizes = slab_sizes,
            fun = apply_fun,
            args = args,
            PACKAGE = "RMINC")
    
    # run the garbage collector...
    gcout <- gc()
    
    result_order <-
      order(results_indexed$inds)
    
    results <- results_indexed$vals[result_order]
    
    if(return_indices)
      results <- list(vals = results, 
                      inds = results_indexed$inds[result_order])
    
    collation_function <- match.fun(collate)
    collated_results <- collation_function(results)
    attr(collated_results, "filenames") <- filenames
    attr(collated_results, "likeVolume") <- filenames[1]
    if(use_mask) attr(collated_results, "mask") <- mask
    
    return(collated_results)
  }

#' @description Can execute any R function at every voxel for a set of MINC files
#' @name mincApply
#' @title Apply arbitrary R function at every voxel
#' @param filenames The MINC files over which to apply the function. Have to be
#' the same sampling.
#' @param function.string The function which to apply. Can only take a single
#' argument, which has to be 'x'. See details and examples.
#' @param mask The filename of a mask - function will only be evaluated
#' inside the mask.
#' @param maskval An integer specifying the value inside the mask where to
#' apply the function. If left blank (the default) then anything
#' above 0.5 will be considered inside the mask. This argument
#' only works for mincApply, not pMincApply.
#' @param reduce Whether to filter masked voxels from the resulting object
#' @details mincApply allows one to execute any R function at every voxel of a
#' set of files. There are two variants: mincApply, which works
#' inside the current R session.
#' Unless the function to be applied takes a single argument a
#' wrapper function has to be written. This is illustrated in the
#' examples; briefly, the wrapper function takes a single argument,
#' called 'x', which will take on the voxel values at every voxel.
#' The function has to be turned into a string; the quote function
#' can be very handy. The output of this function also has to be a
#' simple vector or scalar.
#' Note that interpreted R can be very slow. Mindnumbingly slow. It
#' therefore pays to write optimal functions or, whenever available,
#' use the optimized equivalents. In other words and to give one
#' example, use mincLm rather than applying lm, and if lm has to
#' really be applied, try to use lm.fit rather than plain lm.
#' When using the pbs method, one can also set the options --> TMPDIR,MAX_NODES and WORKDIR
#' @return out The output is a matrix with the same number of rows a the
#' file sizes and the same number of columns as output by the
#' function that was applied. Cast into one of the MINC classes
#' to make writing it out with mincWriteVolume easier.
#' @seealso mincWriteVolume,mincMean,mincSd,mincLm,mincAnova
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' ma <- mincApply(gf$jacobians_fixed_2,quote(mean(x))); 
#' mincWriteVolume(ma, "means.mnc")
#' ### run the whole thing in parallel on the local machine
#' testFunc <- function (x) { return(c(1,2))}
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    workers = 4,
#'                    global = c('gf','testFunc'))
#' mincWriteVolume(pout, "pmincOut.mnc")
#' 
#' ### pbs example 1 (hpf)
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' testFunc <- function (x) { return(c(1,2))}
#' options(TMPDIR="/localhd/$PBS_ID")
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    workers = 4,
#'                    method="pbs",
#'                    global = c('gf','testFunc'))
#'
#' ### pbs example 2 (scinet)
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' testFunc <- function (x) { return(c(1,2))}
#' options(WORKDIR="$SCRATCH")
#' options(MAX_NODES=8)
#' options(TMP_DIR="/tmp")
#' pout <- pMincApply(gf$jacobians_fixed_2, 
#'                    quote(testFunc(x)),
#'                    modules=c("intel","openmpi","R/3.1.1"),
#'                    workers = 4,
#'                    method="pbs",
#'                    global = c('gf','testFunc'))
#'                    
#' #NOTE 1: On SCINET jobs are limited to 32*48 hours of cpu time
#' #NOTE 2: On SCINET the Rmpi library must be compiled for use with R 3.1.1
#' }
#' @export
mincApply <- 
  function(filenames, function.string, mask=NULL, maskval=NULL, reduce=FALSE) {
    
    enoughAvailableFileDescriptors(length(filenames))
    
    if (is.null(maskval)) {
      minmask = 1
      maxmask = 99999999
    }
    else {
      minmask = maskval
      maxmask = maskval
    }
    # Need to get one voxel, x, to test number of values returned from function.string
    x <- mincGetVoxel(filenames, 0,0,0)
    test <- eval(function.string)
    
    results <- .Call("minc2_model",
                     as.character(filenames),
                     matrix(),
                     function.string,
                     NULL,
                     as.double(! is.null(mask)),
                     as.character(mask),
                     as.double(minmask),
                     as.double(maxmask),
                     .GlobalEnv,
                     as.double(length(test)),
                     as.character("eval"), PACKAGE="RMINC");
    
    if (length(test) > 1) {
      if (reduce==TRUE) {
        maskV <- mincGetVolume(mask)
        results <- results[maskV > (minmask-0.5) & maskV < (maxmask+0.5),]
      }
      class(results) <- c("mincMultiDim", "matrix")
    }
    else {
      if (reduce==TRUE) {
        maskV <- mincGetVolume(mask)
        results <- results[maskV > (minmask-0.5) & maskV < (maxmask+0.5)]
      }
      class(results) <- c("mincSingleDim", "numeric")
    }
    attr(results, "likeVolume") <- filenames[1]
    
    # run the garbage collector...
    gcout <- gc()
    
    return(results)
  }

#' Voxel-wise ANOVA
#'  
#' Compute a sequential ANOVA at each voxel
#' @param formula The anova formula. The left-hand term consists of the 
#' MINC filenames over which to compute the models at every voxel.
#' @param data The dataframe which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same 
#' length as the input files. ANOVA will only be computed
#' inside the mask.
#' @details This function computes a sequential ANOVA over a set of files.
#' @return Returns an array with the F-statistic for each model specified by 
#' formula with the following attributes: 
#' \itemize{
#' \item{model}{ design matrix}
#' \item{filenames}{ minc file names input}
#' \item{dimensions}{ dimensions of the statistics matrix}
#' \item{dimnames}{ names of the dimensions for the statistic matrix}
#' \item{stat-type}{ types of statistic used}
#' \item{df}{ degrees of freedom of each statistic}
#' } 
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run an ANOVA at each voxel
#' vs <- mincAnova(jacobians_fixed_2 ~ Sex, gf)
#' }
#' @export
mincAnova <- function(formula, data=NULL, subset=NULL, mask=NULL) {
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  filenames <- as.character(mf[,1])
  mmatrix <- model.matrix(formula, mf)
  
  method <- "anova"
  
  mincFileCheck(filenames)
  
  v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  #l <- lm(formula, mf)
  
  ###   result <- .Call("minc2_model",
  ###                   as.character(filenames),
  ###                   as.matrix(mmatrix),
  ###                   attr(mmatrix, "assign"),
  ###                   as.double(! is.null(mask)),
  ###                   as.character(mask),
  ###                   NULL, NULL,
  ###                   as.character(method), PACKAGE="RMINC")
  result <- .Call("per_voxel_anova",
                  as.character(filenames),
                  as.matrix(mmatrix),
                  attr(mmatrix, "assign"),
                  as.double(! is.null(mask)),
                  as.character(mask), PACKAGE="RMINC")
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- rep("F", ncol(result))
  
  l <- lm.fit(mmatrix, v.firstVoxel)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
    dflist[[i]] <- c(dfs[[i+1]], dfr)
  }
  attr(result, "df") <- dflist
  
  # get the first voxel in order to get the dimension names
  rows <- sub('mmatrix', '',
              rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  
  colnames(result) <- attr(terms(formula), "term.labels")
  class(result) <- c("mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

mincLm_c_wrapper <-
  function(plm, mask, mask_min, mask_max){
    .Call("minc2_model",
          as.character(plm$data.matrix.left),
          plm$data.matrix.right,
          as.matrix(plm$mmatrix),
          NULL,
          as.double(! is.null(mask)),
          as.character(mask),
          as.double(mask_min),
          as.double(mask_max),
          NULL, NULL,
          as.character("lm"), PACKAGE="RMINC")
  }

#' @description Linear Model at Every Voxel
#' @name mincLm
#' @title Linear model at Every Voxel
#' @param formula The linear model formula. The left-hand term consists of the MINC filenames over which to compute the models at every voxel.The RHS of the formula may contain one term with filenames. If so only the + operator may be used, and only two terms may appear on the RHS
#' @param data The data frame which contains the model terms.
#' @param subset Subset definition.
#' @param mask Either a filename or a vector of values of the same length as the input files. The linear model will only be computed
#' inside the mask.
#' @param maskval the value in the mask used to select unmasked voxels, defaults to any positive intensity
#' from 1-99999999 internally expanded to .5 - 99999999.5. If a number is specified voxels with intensities 
#' within 0.5 of the chosen value are considered selected. 
#' @param parallel how many processors to run on (default=single processor).
#' Specified as a two element vector, with the first element corresponding to
#' the type of parallelization, and the second to the number
#' of processors to use. For local running set the first element to "local" or "snowfall"
#' for back-compatibility, anything else will be run with BatchJobs see \link{pMincApply}
#' and \link{configureMincParallel} for details.
#' Leaving this argument NULL runs sequentially.
#' @details This function computes a linear model at every voxel of a set of files. The function is a close cousin to lm, the key difference
#' being that the left-hand side of the formula specification takes a series of filenames for MINC files.
#' @return mincLm returns a mincMultiDim object which contains a series of columns corresponding to the terms in the linear model. The first
#' column is the F-statistic of the significance of the entire volume, the following columns contain the R-Squared term, the marginal t-statistics for each of the terms in the model along with their respective coefficients.
#' @seealso mincWriteVolume,mincFDR,mincMean, mincSd
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # Compute a linear model at each voxel
#' vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
#' }
#' @export
mincLm <- function(formula, data=NULL,subset=NULL , mask=NULL, maskval=NULL, parallel = NULL) {
  
  #INITIALIZATION
  method <- "lm"
  
  # Build model.frame
  m <- m_orig <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  
  if (is.null(maskval)) {
    minmask = 1
    maxmask = 99999999
  }
  else {
    minmask = maskval
    maxmask = maskval
  }
  
  
  # the following returns:
  #
  # list(data.matrix.left = data.matrix.left, 
  #      data.matrix.right = data.matrix.right,
  #      rows = rows,
  #      matrixFound = matrixFound,
  #      mmatrix = mmatrix)
  parseLmOutput <- parseLmFormula(formula,data,mf)
  
  # Call subroutine based on whether matrix was found
  if(parseLmOutput$matrixFound) {
    mincFileCheck(parseLmOutput$data.matrix.left)
    mincFileCheck(parseLmOutput$data.matrix.right)
    enoughAvailableFileDescriptors(length(c(parseLmOutput$data.matrix.left,
                                            parseLmOutput$data.matrix.right)))
  }
  else  {
    parseLmOutput$mmatrix <- model.matrix(formula, mf)	
    parseLmOutput$data.matrix.left <- as.character(mf[,1])
    mincFileCheck(parseLmOutput$data.matrix.left)
    enoughAvailableFileDescriptors(length(parseLmOutput$data.matrix.left))
    parseLmOutput$rows = colnames(parseLmOutput$mmatrix)
  }
  
  # Do Error Checking on mask. The mask may have a different dimension order than the input volumes
  # A call to minc.dimensions.sizes will read the dimension order according to the file order so 
  # it can be used to do the checking.
  if(!is.null(mask)) {
    maskDim = minc.dimensions.sizes(mask)
    volumeDim = minc.dimensions.sizes(parseLmOutput$data.matrix.left)
    # Case 1: Mask has one or more dimensions less than volume
    if(maskDim[1] != volumeDim[1] || maskDim[2] != volumeDim[2] || maskDim[3] != volumeDim[3]) {
      cat(paste("Mask Volume",as.character(maskDim[1]),as.character(maskDim[2]),as.character(maskDim[3]),
                "Input Volume",as.character(volumeDim[1]),as.character(volumeDim[2]),as.character(volumeDim[3]),'\n'))
      stop("Mask Volume input volume dimension mismatch.") }
  }
  
  if(is.null(parallel)){
    result <- mincLm_c_wrapper(parseLmOutput, mask, mask_min = minmask, mask_max = maxmask)
  } else {
    #Setup shared parallel mask and grouping
    n_groups <- as.numeric(parallel[2])
    groups <- seq_len(n_groups)
    new_mask_file <- create_parallel_mask(sample_file = parseLmOutput$data.matrix.left[1]
                                          , mask = mask
                                          , n = n_groups
                                          , tinyMask = tinyMask)
    on.exit(try(unlink(new_mask_file)))
    mask_vol <- as.integer(mincGetVolume(new_mask_file))
    # a vector with two elements: the methods followed by the # of workers
    if (parallel[1] %in% c("local", "snowfall")) {
      result <- 
        parallel::mclapply(groups, function(group){
          mincLm_c_wrapper(parseLmOutput, mask, mask_min = group, mask_max = group)[mask_vol == group, ]
        }, mc.cores = n_groups) %>%
        Reduce(rbind, ., NULL) 
    }
    else {
      reg <- makeRegistry("mincLm_registry")
      on.exit( try(removeRegistry(reg, ask = "no")), add = TRUE)
      
      batchMap(reg, function(group){
        mincLm_c_wrapper(parseLmOutput, mask, mask_min = group, mask_max = group)[mask_vol == group, ]
      }, group = groups)
      
      submitJobs(reg)
      waitForJobs(reg)
      
      result <-
        loadResults(reg, use.names = FALSE) %>%
        Reduce(rbind, ., NULL) 
    } 
    
    result_fleshed_out <- matrix(0
                                 , nrow = prod(minc.dimensions.sizes(new_mask_file))
                                 , ncol = ncol(result))
    
    result_fleshed_out[mask_vol > .5,] <- result 
    result <- result_fleshed_out
  }
  
  attr(result, "likeVolume") <- parseLmOutput$data.matrix.left[1]
  attr(result, "filenames") <- parseLmOutput$data.matrix.left
  attr(result, "model") <- as.matrix(parseLmOutput$mmatrix)
  attr(result, "call") <- m_orig
  attr(result, "data") <- data
  
  
  # the order of return values is:
  #
  # f-statistic
  # r-squared
  # betas
  # t-stats
  #
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))
  
  dflist <- vector("list", (ncol(result)-2)/2 + 1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  # get the first voxel in order to get the dimension names
  #v.firstVoxel <- mincGetVoxel(filenames, 0,0,0)
  #rows <- sub('mmatrix', '',
  #            rownames(summary(lm(v.firstVoxel ~ mmatrix))$coefficients))
  betaNames = paste('beta-',parseLmOutput$rows, sep='')
  tnames = paste('tvalue-',parseLmOutput$rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("mincLm", "mincMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
  
}

#' Minc T-test
#' 
#' Perform an unpaired,unequal variance t-test across a set of minc volumes
#' @param filenames Filenames of the MINC volumes across which to run the t-test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups with which to compare means
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files, with that voxel's t-statistic
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mtt <- mincTtest(gf$jacobians_0.2,gf$Strain)
#' } 
#' @export
mincTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  # the grouping for a t test should only contain 2 groups. Should 
  # also be converted to a factor if it's not.
  if( ! is.factor(grouping) ){
    grouping <- as.factor(grouping)
  }
  if(length(levels(grouping)) != 2 ){
    cat("\nThe groups (levels) in your data are:\n")
    cat(levels(grouping), "\n\n")
    stop("mincTtest can only be performed on data with 2 groups")
  }
  result <- mincSummary(filenames, grouping, mask, method="t-test", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("t") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  rttest <- t.test(vox~grouping,data=gf,paired=FALSE)
  attr(result, "df") <- rttest$parameter
  colnames(result) <- c("T-statistic")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)
}

#' @title Minc Paired T Test
#' @description Perform a paired t-test across a set of minc volumes
#' @param filenames Filenames of the MINC volumes across which to run the t-test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups with which to compare means. The two groups
#' must be the same length.
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#'          elements as there are voxels in the input files, with that voxel's t-statistic
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' gf = gf[1:20,]
#' mptt <- mincPairedTtest(gf$jacobians_0.2,gf$Strain)
#' }
#' @export
mincPairedTtest <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  # here, similarly to the t-test, there should be 2 groups in the data. However
  # since this is a paired t-test (repeated measures, and it assumes element 1 from
  # group 1 is paired with element 1 from group 2 etc.), the groups need to have 
  # the same length.
  if( ! is.factor(grouping) ){
    grouping <- as.factor(grouping)
  }
  if(length(levels(grouping)) != 2 ){
    cat("\nThe groups (levels) in your data are:\n")
    cat(levels(grouping), "\n\n")
    stop("mincPairedTtest can only be performed on data with 2 groups")
  }
  lenght_group_1 <- sum(grouping == levels(grouping)[1])
  lenght_group_2 <- sum(grouping == levels(grouping)[2])
  if(lenght_group_1 != lenght_group_2) {
    cat("\nGroup 1 (",levels(grouping)[1], ") has ", lenght_group_1," elements.\n")
    cat("Group 2 (",levels(grouping)[2], ") has ", lenght_group_2," elements.\n\n")
    stop("The groups do not have the same length. mincPairedTtest expects 2 groups with equal length (paired data)")
  }
  result <- mincSummary(filenames, grouping, mask, method="paired-t-test", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("t") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  rttest <- t.test(vox~grouping,data=gf,paired=TRUE)
  attr(result, "df") <- rttest$parameter
  colnames(result) <- c("T-statistic")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)
}

#' @title Minc Wilcoxon
#' @description Perform a Mann-Whitney U test between a set of minc volumes.
#' @param filenames Filenames of the MINC volumes across which to run the test
#' @param grouping  Contains same number of elements as
#' filenames; must contain exactly two groups.
#' @param mask A mask specifying which voxels are to be included in the test
#' @param maskval The value with which to mask the data (data will masked +/- 0.5 around this value
#' @return  The output will be a single vector containing as many
#' elements as there are voxels in the input files, with that voxel's U value (lower one).
#' The number of observations in each sample is also saved as an attribute(m and n) so the result
#' can be passed into mincFDR.
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
#' mw <- mincWilcoxon(gf$jacobians_0.2,gf $Strain)
#' }
#' @export
mincWilcoxon <- function(filenames, grouping, mask=NULL, maskval=NULL) {
  result <- mincSummary(filenames, grouping, mask, method="wilcoxon", maskval=maskval)
  result <- as.matrix(result);
  attr(result, "likeVolume") <- filenames[1]
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- c("u") 
  gf <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
  gf$grouping <- grouping
  gf$vox <- mincGetVoxel(filenames, 0, 0, 0)
  a = levels(gf$grouping)
  attr(result, "m") <- nrow(gf[gf$grouping == a[1],])
  attr(result, "n") <- nrow(gf[gf$grouping == a[2],])
  colnames(result) <- c("Mann-Whitney")
  class(result) <- c("mincMultiDim", "matrix")
  return(result)
  return(result)
}


#' Threshold Free Cluster Enhancement
#' 
#' Perform threshold free cluster enhancement as described in 
#' Smith and Nichols (2008). Cluster-like structures are enhanced
#' to allow a hybrid cluster/voxel analysis to be performed.
#' @param volume Either a character vector with a single filename, 
#' a \code{mincSingleDim}, or a \code{mincMultiDim} object
#' @param d The discretization step-size for approximating the threshold integral (default .1)
#' @param E The exponent by which to raise the extent statistic (default .5)
#' @param H The exponent by which to raise the height (default 2)
#' @param side Whether to consider positive and negative statistics or both (default both)
#' @param output_file A filename for the enhanced volume.
#' @param keep whether to keep the resulting volume, defaults to whether or not you specified
#' a name for the output volume.
#' @param ... additional arguments for methods
#' @return A minc volume with the results of the enhancement
#' @export  
mincTFCE <-
  function(volume, d = 0.1, E = .5, H = 2.0
           , side = c("both", "positive", "negative")
           , output_file = NULL
           , keep = is.null(output_file)
           , ...) {
    
    UseMethod("mincTFCE")
  }

#' @method mincTFCE mincSingleDim
#' @export 
mincTFCE.mincSingleDim <-
  function(volume, d = 0.1, E = .5, H = 2.0
           , side = c("both", "positive", "negative")
           , output_file = NULL
           , keep = is.null(output_file)
           , ...){
    
    side <- match.arg(side)
    in_file <- tempfile("minc_tfce_in", fileext = ".mnc")
    on.exit(unlink(in_file))
    
    if(is.null(output_file))
      output_file <- tempfile("minc_tfce_out", fileext = ".mnc")
    
    if(!keep)  on.exit(unlink(output_file), add = TRUE)

    capture.output(mincWriteVolume(volume, in_file))
    
    side_flag <-
      `if`(side == "both"
           , "pos-and-neg"
           , `if`(side == "positive"
                  , "pos-only"
                  , "neg-only"))
    
    system(
      sprintf("TFCE -d %s -E %s -H %s --%s %s %s"
              , d, E, H
              , side_flag
              , in_file
              , output_file)
      , intern = TRUE)

      
    out_vol <- mincGetVolume(output_file)
    likeVolume(out_vol) <- likeVolume(volume)
    
    out_vol
  }

mincTFCE.mincLm <- 
  function(model, tails = c("both", "positive", "negative")
           , d = 0.1, E = .5, H = 2.0
           , side = c("both", "positive", "negative")
           , output_file = NULL
           , ...){
    m_call <- attr(model, "call")
    m_data <- attr(model, "model")
    
    random_data <- m_data[sample(seq_len(nrow(data))),]
    
    m_call["data"] <- random_data
    
    lm_res <- eval.parent(m_call)
    
    
  }


