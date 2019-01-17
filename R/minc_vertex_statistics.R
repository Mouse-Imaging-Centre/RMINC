#' Create descriptive statistics across a series of vertex files
#' 
#' This function is used to compute the mean, standard deviation, sum, or
#' variance of every vertex in a set of vertex files.
#' 
#' 
#' @param filenames Filenames of the vertex volumes across which to create the
#' descriptive statistic.
#' @return \item{out}{The output will be a single vector containing as many
#' elements as there are vertices in the input files.}
#' @seealso vertexLm
#' @examples
#' 
#' \dontrun{
#' # read the text file describing the dataset
#' gf <- read.csv("control-file.csv")
#' # compute the mean at every voxel of all files.
#' means <- vertexMean(gf$filenames)
#' }
#' @name vertexSummaries
NULL

#' @describeIn vertexSummaries mean
#' @export
vertexMean <- function(filenames, column=1) 
{
  vertexData = vertexTable(filenames, column=column)
  return(rowMeans(vertexData))
  
} 

#' @describeIn vertexSummaries sum
#' @export
vertexSum <- function(filenames, column=1) 
{
  vertexData = vertexTable(filenames,column=column)
  return(rowSums(vertexData))
  
} 

#' @describeIn vertexSummaries var
#' @export
vertexVar <- function(filenames, column=1) 
{
  vertexData = vertexTable(filenames,column=column)
  return(apply(vertexData,1,var))
  
} 

#' @describeIn vertexSummaries standard deviation
#' @export
vertexSd<- function(filenames,column=1) 
{
  vertexData = vertexTable(filenames,column=column)
  return(apply(vertexData,1,sd))
  
}

### Helper function for applying over rows of a potentially 
### masked matrix potentially in parallel
matrixApply <- function(mat, fun, ..., mask = NULL, parallel = NULL
                      , collate = simplify_masked
                      , resources = list() 
                      , conf_file = getOption("RMINC_BATCH_CONF")){
  
  if(!is.null(mask)){
    if(length(mask) == 1 && is.character(mask))
      mask <- readLines(mask)
    
    mask_lgl <- mask > .5
    mat <- mat[mask_lgl,]
  }
  
  fun <- match.fun(fun)
  
  apply_fun <- function(sub_matrix){
    if(!is.matrix(sub_matrix))
      sub_matrix <- matrix(sub_matrix, nrow = 1)
    
    lapply(seq_len(nrow(sub_matrix))
           , function(i) fun(sub_matrix[i,], ...))
  }
  
  if(is.null(parallel)){
    results <- apply_fun(mat)
  } else {
    n_groups <- as.numeric(parallel[2])
    groups <- split(seq_len(nrow(mat)), groupingVector(nrow(mat), n_groups))
    
    if(parallel[1] == "local") {
      results <- 
        failing_mclapply(groups, function(group){
          apply_fun(mat[group,])
        }, mc.cores = n_groups) %>%
        Reduce(c, ., NULL)
    } else {
      reg <- qMincRegistry(new_file("matrixApply_registry"), conf_file = conf_file)
      on.exit( tenacious_remove_registry(reg) )
      
      ids <- batchMap(reg = reg, function(group){
        apply_fun(mat[group,])
      }, group = groups)
      
      submitJobs(ids, reg = reg, resources = resources)
      waitForJobs(reg = reg)
      
      results <-
        reduceResults(c, init = NULL, reg = reg)
    }
  }
  
  # Flesh out the object if a mask was used
  if(!is.null(mask)){
    mask_val <- getOption("RMINC_MASKED_VALUE")
    results_expanded <- replicate(length(mask_lgl), getOption("RMINC_MASKED_VALUE"), simplify = FALSE)
    results_expanded[mask_lgl] <- results
    results <- results_expanded
  }
  
  collate_fun <- match.fun(collate)
  results <- collate(results)
  
  results
}

#' @description This function is used to compute an arbitrary function of every region in a set of vertex files.
#' @name vertexApply
#' @title Apply function over vertex Files
#' @param filenames vertex file names
#' @param fun A function to be applied to each vertex
#' @param ... additional arguments to \code{fun}
#' @param mask A vector of filename indicating a vertex mask (\code{fun} is applied to all vertices
#' where mask is greater than .5)
#' @param parallel A two component vector indicating how to parallelize the computation. If the 
#' first element is "local" the computation will be run via the parallel package, otherwise it will
#' be computed using batchtools, see \link{pMincApply} for details. The element should be numeric
#' indicating the number of jobs to split the computation into.
#' @param collate A function to reduce the (potentially masked) list of results into a nice
#' structure. Defaults to \link{simplify_masked}
#' @param transpose Whether to alternatively transpose the vertex matrix and apply a function to
#' each subject
#' @return  The a matrix with a row of results for each vertex
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- vertexApply(gf$CIVETFILES$nativeRMStlink20mmleft, mean)
#' }
#' @export
vertexApply <- function(filenames, fun, ..., mask = NULL, parallel = NULL
                      , collate = simplify_masked, transpose = FALSE, column=1) 
{
  # Load the data
  vertexData <- vertexTable(filenames, column=column)
  if(transpose)
    vertexData <- t(vertexData)
  
  results <- matrixApply(vertexData, fun, ..., mask = mask, parallel = parallel, collate = collate)
  attr(results, "likeFile") <- filenames[1]
  
  results
}


#' Performs ANOVA on each vertex point specified 
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param subset rows to be used, by default all are used
#' @return Returns an array with the F-statistic for each model specified by formula with the following attributes: 
#' \itemize{
#' \item{model}{ design matrix}
#' \item{filenames}{ minc file names input}
#' \item{dimensions}{ dimensions of the statistics matrix}
#' \item{dimnames}{ names of the dimensions for the statistic matrix}
#' \item{stat-type}{ types of statistic used}
#' \item{df}{ degrees of freedom of each statistic}
#' } 
#' @seealso mincAnova,anatAnova 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' result = vertexAnova(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf)
#' }
#' @export 
vertexAnova <- function(formula, data, subset=NULL, column=1) {
  # Create Model
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mmatrix <- model.matrix(formula, mf)
  
  # Load Vertex Data from Files
  filenames <- as.character(mf[,1])
  data.matrix <- vertexTable(filenames, column=column)
  result <- .Call("vertex_anova_loop", data.matrix, mmatrix,attr(mmatrix, "assign"), PACKAGE="RMINC");
  
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "filenames") <- filenames
  attr(result, "stat-type") <- rep("F", ncol(result))
  
  # Use the first voxel in order to get the dimension names
  v.firstVoxel <- data.matrix[1,]
  
  # Get the degrees of freedom
  l <- lm.fit(mmatrix, v.firstVoxel)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
    dflist[[i]] <- c(dfs[[i + 1]], dfr)
  }
  attr(result, "df") <- dflist
  
  colnames(result) <- attr(terms(formula), "term.labels")
  
  class(result) <- c("vertexMultiDim", "matrix")
  return(result)
}


#' Calculates statistics and coefficients for linear model of specified vertex files
#' @param formula a model formula. The RHS of the formula may contain one term with filenames. If
#' so only the + operator may be used, and only two terms may appear on the RHS
#' @param data a data.frame containing variables in formula 
#' @param subset rows to be used, by default all are used
#' @return Returns an object containing the R-Squared value,beta coefficients, F 
#' and t statistcs that can be passed directly into vertexFDR.
#' @seealso mincLm,anatLm,vertexFDR 
#' @examples
#' \dontrun{ 
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' result = vertexLm(CIVETFILES$nativeRMStlink20mmleft~Primary.Diagnosis,gf) 
#' }
#' @export
vertexLm <- function(formula, data, subset=NULL, column=1 ) {
  # Build model.frame
  m <- m_orig <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)

  mf <- mf[c(1, m)]
  
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  
  mf <- eval(mf, parent.frame())

  mincFileCheck(as.character(mf[,1]))

  if(length(grep("\\$",formula[[3]])) > 0) {
    stop("$ Not Permitted in Formula")  
  }
  
  # the following returns:
  #
  # list(data.matrix.left = data.matrix.left, 
  #      data.matrix.right = data.matrix.right,
  #      rows = rows,
  #      matrixFound = matrixFound,
  #      mmatrix = mmatrix)
  parseLmOutput <- parseLmFormula(formula,data,mf)

  if(parseLmOutput$matrixFound) {

    parseLmOutput$data.matrix.left  <- vertexTable(parseLmOutput$data.matrix.left, column=column)
    parseLmOutput$data.matrix.right <- vertexTable(parseLmOutput$data.matrix.right, column=column)
  }
  else {

    filenames <- as.character(mf[,1])
    parseLmOutput$mmatrix <- model.matrix(formula, mf)	
    parseLmOutput$data.matrix.left <- vertexTable(filenames, column=column)
    parseLmOutput$rows = colnames(parseLmOutput$mmatrix)
  }
  result <- .Call("vertex_lm_loop",
                  parseLmOutput$data.matrix.left,
                  parseLmOutput$data.matrix.right,
                  parseLmOutput$mmatrix,
                  PACKAGE="RMINC") 

  attr(result, "likeVolume") <- as.character(mf[,1])[1]
  attr(result, "model")      <- as.matrix(parseLmOutput$mmatrix)
  attr(result, "filenames")  <- as.character(mf[,1])
  attr(result, "stat-type")  <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2), "logLik")
  attr(result, "data")       <- data 
  attr(result, "call")       <- m_orig

  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))
  
  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2+1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  betaNames = paste('beta-', parseLmOutput$rows, sep='')
  tnames = paste('tvalue-', parseLmOutput$rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames, "logLik")
  class(result) <- c("vertexLm", "vertexMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  return(result)
}

#' Vertex Mixed Effects Models
#' 
#' Perform linear mixed effects model fitting for vertex data.
#' vertexLmer should be used the same way as a straight lmer call, except
#' that the left hand side of the equation contains vertex filenames rather than
#' an actual response variable.
#' 
#' @inheritParams mincLmer
#' @details \code{vertexLmer}, like its relative \link{mincLmer} provides an interface to running 
#' linear mixed effects models at every vertex. Unlike standard linear models testing hypotheses 
#' in linear mixed effects models is more difficult, since the denominator degrees of freedom are 
#' more difficult to  determine. RMINC provides estimating degrees of freedom using the
#' \code{\link{vertexLmerEstimateDF}} function. For the most likely models - longitudinal
#' models with a separate intercept or separate intercept and slope per subject - this
#' approximation is likely correct. Be careful in using this approximation if
#' using more complicated random effects structures.
#'
#' @seealso \code{\link{lmer}} for description of lmer and lmer formulas; \code{\link{mincLm}}
#' @export
vertexLmer <-
  function(formula, data, mask=NULL, parallel=NULL,
           REML=TRUE, control=lmerControl(), start=NULL,
           verbose=0L, safely = FALSE, summary_type = "fixef") {

    mc <- mcout <- match.call()
    mc[[1]] <- quote(lme4::lFormula)

    # remove lme4 unknown arguments, since lmer does not know about them and keeping them
    # generates obscure warning messages
    mc <- mc[!names(mc) %in% c("mask", "parallel", "safely", "summary_type")]

    lmod <- eval(mc, parent.frame(1L))
    mincFileCheck(lmod$fr[,1])
    if(!is.null(mask) & (!(is.character(mask) & length(mask) == 1)) & !is.numeric(mask))
      stop("A numeric vector or file name must be supplied for the mask argument")

    # code ripped from lme4:::mkLmerDevFun
    rho <- new.env(parent = parent.env(environment()))
    rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta",
                                                    "Lambdat", "Lind")],
                                      n = nrow(lmod$X), list(X = lmod$X)))
    REMLpass <- if (REML)
      ncol(lmod$X)
    else 0L


    mincLmerList <- list(lmod, mcout, control, start, verbose, rho, REMLpass)

    summary_fun <- summary_type
    if(is.character(summary_type) && length(summary_type) == 1)
       summary_fun <- switch(summary_type
                           , fixef = fixef_summary
                           , ranef = ranef_summary
                           , both = effect_summary
                           , anova = anova_summary
                           , stop("invalid summary type specified"))

    if(!is.function(summary_fun))
      stop("summary_type must be a string specifying a summary, or a function")
    

    mincLmerOptimizeAndExtractSafely <-
      function(x, mincLmerList, summary_fun){
        tryCatch(mincLmerOptimizeAndExtract(x, mincLmerList, summary_fun),
                 error = function(e){warning(e); return(NA)})
      }

    optimizer_fun <-
      `if`(safely, mincLmerOptimizeAndExtractSafely, mincLmerOptimizeAndExtract)

    out <- 
      vertexApply(lmod$fr[,1]
                  , optimizer_fun
                  , mincLmerList = mincLmerList
                  , mask = mask
                  , summary_fun = summary_fun
                  , parallel = parallel)

    ## Result post processing
    out[is.infinite(out)] <- 0            #zero out infinite values produced by vcov

    # generate some random numbers for a single fit in order to extract some extra info
    mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])), mincLmerList)

    res_cols <- colnames(out)
    attr(out, "stat-type") <- ## Handle all possible output types
      check_stat_type(res_cols, summary_type)

    # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
    attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
    attr(out, "REML") <- REML
    attr(out, "mask") <- mask
    attr(out, "data") <- data
    attr(out, "mincLmerList") <- mincLmerList
    class(out) <- c("vertexLmer", "mincLmer", "mincMultiDim", "matrix")

    return(out)
  }

#' Estimate the degrees of freedom for parameters in a vertexLmer model
#'
#' There is much uncertainty in how to compute p-values for mixed-effects
#' statistics, related to the correct calculation of the degrees of freedom
#' of the model (see here \url{http://glmm.wikidot.com/faq#df}). mincLmer by
#' default does not return the degrees of freedom as part of its model, instead
#' requiring an explicit call to a separate function (such as this one).
#' The implementation here is the Satterthwaite approximation. This approximation
#' is computed from the data, to avoid the significant run-time requirement of computing
#' it separate for every vertex, here it is only computed on a small number of vertices
#' within the mask and the median DF returned for every variable.
#' 
#' @param model the output of mincLmer
#'
#' @return the same mincLmer model, now with degrees of freedom set
#'
#' @seealso \code{\link{mincLmer}} for mixed effects modelling, \code{\link{mincFDR}}
#' for multiple comparisons corrections.
#'
#' @examples
#' \dontrun{
#' vs <- mincLmer(filenames ~ age + sex + (age|id), data=gf, mask="mask.mnc")
#' vs <- mincLmerEstimateDF(vs)
#' qvals <- mincFDR(vs, mask=attr(vs, "mask"))
#' qvals
#' }
#' @export
vertexLmerEstimateDF <-
  function(model, column=1){
    # set the DF based on the Satterthwaite approximation
    mincLmerList <- attr(model, "mincLmerList")
    mask <- attr(model, "mask")
    
    if(is.null(mask)){
      mask <- rep(1, nrow(model))
    } else if(is.character(mask) && length(mask) == 1) {
      mask <- as.numeric(readLines(mask))
    } else if(! (is.numeric(mask) || is.logical(mask)) || length(mask) != nrow(model)){
      stop("There is a problem with your mask, please check that it fits your results")
    }
    
    initial_frame <- #unpack the lmerList object to get the raw data
      attr(model, "mincLmerList")[[1]]$fr
    
    vertex_data <- vertexTable(initial_frame[,1], column=column)
    
    # estimated DF depends on the input data. Rather than estimate separately at every structure,
    # instead select a small number of structures and estimate DF for those structures, then keep the
    # min
    nvertices <- min(50, sum(mask > .5))
    rvertices <- sample(which(mask > .5), nvertices)
    dfs <- matrix(nrow = nvertices, ncol = sum(attr(model, "stat-type") %in% "tlmer"))
    original_data <- attr(model, "data")
    
    ## It seems LmerTest cannot compute the deviance function for mincLmers
    ## in the current version, instead extract the model components from
    ## the mincLmerList and re-fit the lmers directly at each structure,
    ## Slower but yeilds the correct result
    lmod <- mincLmerList[[1]]
    LHS <- as.character(lmod$formula[[2]])
    form <- update(lmod$formula, RMINC_DUMMY_LHS ~ .)

    filenames <- unique(mincLmerList[[1]]$fr[,1])
    row_file_match <- match(original_data[[LHS]], filenames)

    for (i in 1:nvertices) {
      vertex_vals <- vertex_data[rvertices[i],]
      original_data$RMINC_DUMMY_LHS <- vertex_vals[row_file_match]

      model_env <- list2env(original_data)
      environment(form) <- model_env

      mmod <-
        lmerTest::lmer(form, REML = lmod$REML,
                       start = mincLmerList[[4]], control = mincLmerList[[3]],
                       verbose = mincLmerList[[5]])

      dfs[i,] <- 
        suppressMessages(
          tryCatch(summary(mmod)$coefficients[,"df"]
                 , error = function(e){ 
                   warning("Unable to estimate DFs for vertex "
                         , rvertices[i]
                         , call. = FALSE)
                   NA}))
    }
    
    df <- apply(dfs, 2, median, na.rm = TRUE)
    cat("Mean df: ", apply(dfs, 2, mean, na.rm = TRUE), "\n")
    cat("Median df: ", apply(dfs, 2, median, na.rm = TRUE), "\n")
    cat("Min df: ", apply(dfs, 2, min, na.rm = TRUE), "\n")
    cat("Max df: ", apply(dfs, 2, max, na.rm = TRUE), "\n")
    cat("Sd df: ", apply(dfs, 2, sd, na.rm = TRUE), "\n")
    
    attr(model, "df") <- df
    
    return(model)
  }

#' Threshold Free Cluster Enhancement
#' 
#' Perform threshold free cluster enhancement as described in 
#' Smith and Nichols (2008). Cluster-like structures are enhanced
#' to allow a hybrid cluster/voxel analysis to be performed.
#' @param x A numeric vector, a filepath to a set of values,
#' or a \code{matrix} object, or \code{vertexLm} object.
#' @param surface Either a mesh object corresponding to the surface,
#' an igraph graph object of surface created by \link{obj_to_graph}, or an adjacency list (see details). 
#' For the \code{matrix} and \link{vertexLm} cases, either a single surface object may be passed and
#' used for each individual, or a vector of file names
#' @param nsteps The number of steps to discretize the TFCE computation over
#' @inheritParams mincTFCE
#' @param weights A weighting vector assigning area to vertices. The default varies by 
#' \code{surface} object type. If \code{surface} is an \code{bic_obj} the area of each
#' vertex is equal to 1/3 the sum of the areas of the triangles it is a member of. If \code{surface}
#' is an \code{igraph} or adjacency list the default is a vector of ones.
#' @details Passing an adjacency list will save some compute time
#' but is not recommended for general use. If an adjacency list is passed should index starting from 0 for compatibility with c++
#' code. Adjacency lists of this kind can be generated from graphs with \code{lapply(as_adj_list(graph), `-`, 1)}
#' using the \link[igraph]{as_adj_list} from the igraph library. 
#' @return The behaviour of \code{vertexTFCE} is to perform cluster free enhancement on a object,
#' in the single dimensional case, a string denoting a vertex file or a numeric vector
#' it returns a numeric vector with the result. In the matrix case each column is cluster enhanced and 
#' recomposed into a matrix.
#' In the \link{vertexLm} case a randomization test is performed on each t-tstatistic column. The results is
#' a list with 3 elements \itemize{
#' \item{TFCE: A matrix of the tvalue columns after randomization}
#' \item{randomization_dist: An RxT matrix where R is the number of randomizations and T is the number of t-statistics,
#' elements are the largest value obtained by the randomized TFCE}
#' \item{args: Arguments passed to the internal randomzations and TFCE code}
#' }
#' @export  
vertexTFCE <-
  function(x, ...) {
    
    UseMethod("vertexTFCE")
  }

#' @describeIn vertexTFCE numeric
#' @export 
vertexTFCE.numeric <-
  function(x, surface, E = .5, H = 2.0
           , nsteps = 100
           , side = c("both", "positive", "negative")
           , weights = NULL
           , ...){
     side <- match.arg(side)
     
     if(is.null(weights)){
       if(inherits(surface, "bic_obj")){
         weights <- mesh_area(surface$vertex_matrix, surface$triangle_matrix)
       } else {
         weights <- rep(1, length(x))
       }
     } else {
       if(length(weights) == 1) weights <- rep(weights, length(x))
       if(length(weights) != length(x)) stop("Please supply 1 or enough weights for each element of x")
     }
     
     if(inherits(surface, "character")){
       if(length(surface) != 1) stop("Only one file can be passed for the surface")
       surface <- read_obj(surface)
     }
     
     if(inherits(surface, "bic_obj")){
       surface <- obj_to_graph(surface)
     }
     
     if(inherits(surface, "igraph")){
       igraph_return_vs_old <- igraph::igraph_opt("return.vs.es")
       on.exit(igraph::igraph_options(return.vs.es = igraph_return_vs_old))
       igraph::igraph_options(return.vs.es = FALSE)
       adjacencies <- lapply(igraph::as_adj_list(surface), `-`, 1)
     } else {
       adjacencies <- surface
     }
       
     
     if(length(weights) == 1)
       weights <- rep(weights, length(adjacencies))
     
     
     tfce <-
       switch(side
              , "positive" = graph_tfce_wqu(x, adjacencies, E, H, nsteps, weights)
              , "negative" = -graph_tfce_wqu(-x, adjacencies, E, H, nsteps, weights)
              , "both" = graph_tfce(x, adjacencies, E, H, nsteps, weights))
     
     return(tfce)
  }

#'@describeIn vertexTFCE matrix
#'@export
vertexTFCE.matrix <-
  function(x, surface, E = .5, H = 2.0
           , nsteps = 100
           , side = c("both", "positive", "negative")
           , weights = NULL
           , ...){
    
    side <- match.arg(side)
    length_err <- "You must pass either a single surface or one per column of the input matrix"
    
    ## Irritating case analysis - handle lists of object,graphs, adjacency lists, or character vectors,
    ## or singletons of each of those.
    if(is.list(surface) && inherits(surface[[1]], c("bic_obj", "igraph"))){
      ## Case: list of objects or graphs
      if(length(surface) != ncol(x))  stop(length_err)
      surface_inds <- seq_along(surface)
      
    } else if(is.list(surface) && 
              !inherits(surface, c("bic_obj", "igraph")) && 
              !inherits(surface[[1]], c("bic_obj", "igraph", "numeric", "integer"))){
      ## Case list of adjacency lists
      if(length(surface) != ncol(x)) stop(length_err)
      surface_inds <- seq_along(surface)
      
    } else if(is.character(surface)){
      ## Case: vector or string of filenames
      if(length(surface) == 1){ 
        surface_inds <- rep(1, ncol(x))
      } else if(length(surface) == ncol(x)){
        surface_inds <- seq_along(surface)
      } else {
        stop(length_err)
      }
      
    } else {
      ## Case: single surface, graph, or adjacency list (convert to an adjacency list to save computation)
      if(is.list(surface) && !inherits(surface, c("bic_obj", "igraph"))) 
        surface <- list(surface)
      
      if(inherits(surface, "bic_obj")) 
        surface <- obj_to_graph(surface)
      
      if(inherits(surface, "igraph")){
        igraph_return_vs_old <- igraph::igraph_opt("return.vs.es")
        on.exit(igraph::igraph_options(return.vs.es = igraph_return_vs_old))
        igraph::igraph_options(return.vs.es = FALSE)
        adjacencies <- lapply(igraph::as_adj_list(surface), `-`, 1)
        surface <- list(surface)
      }
      
      surface_inds <- rep(1, ncol(x))
    }
    
    ## Finally, the TFCE part, passes off most of the work to vertexTFCE.numeric
    tf_res <- 
      mapply(function(col_ind, surface_ind){
        vertexTFCE.numeric(x[,col_ind], surface[[surface_ind]], E = E, H = H, nsteps = nsteps, weights = weights)
      }, seq_len(ncol(x)), surface_inds)
    
    colnames(tf_res) <- colnames(x)
    tf_res
  }

#' @export
getCall.vertexLm <-
  function(x, ...) attributes(x)$call

#'@describeIn vertexTFCE vertexLm
#'@export
vertexTFCE.vertexLm <- 
  function(x
           , surface
           , R = 500
           , alternative = c("two.sided", "greater")
           , E = .5, H = 2.0, nsteps = 100
           , weights = NULL
           , side = c("both", "positive", "negative")
           , replace = FALSE
           , parallel = NULL
           , ...){
    
    lmod          <- x
    alternative   <- match.arg(alternative)
    side          <- match.arg(side)
    columns       <- grep("tvalue-", colnames(lmod))
    lmod_call     <- attr(lmod, "call")
    like_vol      <- likeVolume(lmod) 
    
    original_tfce <-
      vertexTFCE.matrix(lmod[,columns], surface = surface, weights = weights
                        , nsteps = nsteps, E = E, H = H, side = side) 
    
    randomization_dist <-
      mincRandomize_core(lmod, R = R, replace = replace, parallel = parallel, columns = columns
                         , alternative = alternative
                         , post_proc = vertexTFCE.matrix
                         , surface = surface
                         , weights = weights
                         , side = side, E = E, H = H, nsteps = nsteps)
    
    output <- list(tfce = original_tfce, randomization_dist = randomization_dist
                   , args = list(call = lmod_call,
                                 side = side
                                 , alternative = alternative))
    class(output) <- c("vertexTFCE_randomization", "vertex_randomization")
    output
  }

#'@describeIn vertexTFCE character
#'@export
vertexTFCE.character <-
  function(x, surface, E = .5, H = 2.0
           , nsteps = 100
           , side = c("both", "positive", "negative")
           , weights = NULL
           , column=1
           , ...){
    if(length(x) == 1){
      x <- as.numeric(readLines(x))
      vertexTFCE.numeric(x, surface, E = E, H = H, nsteps = nsteps, side = side, weights = weights, ...)
    } else {
      x <- vertexTable(x,column=column)
      vertexTFCE.matrix(x, surface, E = E, H = H, nsteps = nsteps, side = side, weights = weights, ...)
    }
  }

vertexRandomize.vertexLm <-
  function(x
           , R = 500
           , alternative = c("two.sided", "greater")
           , replace = FALSE
           , parallel = NULL
           , columns = grep("tvalue-", colnames(x))){
    lmod          <- x
    alternative   <- match.arg(alternative)
    original_stats <- lmod[,columns]
    lmod_call   <- attr(lmod, "call")
    
    randomization_dist <-
      mincRandomize_core(lmod, R = R, replace = replace, parallel = parallel, columns = columns
                         , alternative = alternative)
    
    output <- list(stats = original_stats, randomization_dist = randomization_dist
                   , args = c(call = lmod_call, alternative = alternative))
    class(output) <- c("vertexLm_randomization", "vertex_randomization")
    
    output
  
  }


# vertexLmer <-
#   function(formula, data, mask=NULL, parallel=NULL,
#            REML=TRUE, control=lmerControl(), start=NULL, 
#            verbose=0L, temp_dir = getwd(), safely = FALSE, 
#            cleanup = TRUE) {
#     
#     mc <- mcout <- match.call()
#     mc[[1]] <- quote(lme4::lFormula)
#     
#     # remove lme4 unknown arguments, since lmer does not know about them and keeping them
#     # generates obscure warning messages
#     mc <- mc[!names(mc) %in% c("mask", "parallel", "temp_dir", "safely", "cleanup")]
#     
#     lmod <- eval(mc, parent.frame(1L))
#     
#     # code ripped from lme4:::mkLmerDevFun
#     rho <- new.env(parent = parent.env(environment()))
#     rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta", 
#                                                     "Lambdat", "Lind")],
#                                       n = nrow(lmod$X), list(X = lmod$X)))
#     REMLpass <- if (REML) 
#       ncol(lmod$X)
#     else 0L
#     
#     
#     mincLmerList <- list(lmod, mcout, control, start, verbose, rho, REMLpass)
#     
#     
#     mincLmerOptimizeAndExtractSafely <-
#       function(x, mincLmerList){
#         tryCatch(mincLmerOptimizeAndExtract(x, mincLmerList),
#                  error = function(e){warning(e); return(NA)})
#       }
#     
#     optimizer_fun <- 
#       `if`(safely, mincLmerOptimizeAndExtractSafely, mincLmerOptimizeAndExtract)
#     
#     if (!is.null(parallel)) {
#       # a vector with two elements: the methods followed by the # of workers
#       if (parallel[1] %in% c("local", "snowfall")) {
#         out <- mcMincApply(lmod$fr[,1],
#                            optimizer_fun,
#                            mincLmerList = mincLmerList,
#                            filter_masked = TRUE,
#                            mask = mask,
#                            cores = as.numeric(parallel[2]),
#                            slab_sizes = slab_dims,
#                            cleanup = cleanup)
#       }
#       else if(parallel[1] %in% c("sge", "pbs")){
#         reg <- qMincRegistry("qVertexLmer_registry",
#                              parallel_method = parallel[1],
#                              temp_dir = temp_dir,
#                              cores = 1)
#         
#       } else {
#         stop("Error: unknown parallelization method")
#       }
#     }
#     else {
#       out <- mincApplyRCPP(lmod$fr[,1], # assumes that the formula was e.g. filenames ~ effects
#                            optimizer_fun,
#                            mincLmerList = mincLmerList,
#                            mask = mask,
#                            slab_sizes = slab_dims)
#     }
#     
#     ## Result post processing
#     out[is.infinite(out)] <- 0            #zero out infinite values produced by vcov
#     
#     termnames <- colnames(lmod$X)
#     betaNames <- paste("beta-", termnames, sep="")
#     tnames <- paste("tvalue-", termnames, sep="")
#     colnames(out) <- c(betaNames, tnames, "logLik", "converged")
#     
#     # generate some random numbers for a single fit in order to extract some extra info
#     mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])), mincLmerList)
#     
#     attr(out, "stat-type") <- c(rep("beta", length(betaNames)), rep("tlmer", length(tnames)),
#                                 "logLik", "converged")
#     # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
#     attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
#     attr(out, "REML") <- REML
#     attr(out, "mask") <- mask
#     attr(out, "mincLmerList") <- mincLmerList
#     class(out) <- c("mincLmer", "mincMultiDim", "matrix")
#     
#     return(out)
#     
#     # out <- t(apply(anat, 2, mincLmerOptimizeAndExtract, mincLmerList = mincLmerList))
#     # 
#     # out[is.infinite(out)] <- 0            #zero out infinite values produced by vcov
#     # 
#     # termnames <- colnames(lmod$X)
#     # betaNames <- paste("beta-", termnames, sep="")
#     # tnames <- paste("tvalue-", termnames, sep="")
#     # colnames(out) <- c(betaNames, tnames, "logLik", "converged")
#     # 
#     # # generate some random numbers for a single fit in order to extract some extra info
#     # mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])), mincLmerList)
#     # 
#     # attr(out, "stat-type") <- c(rep("beta", length(betaNames)), rep("tlmer", length(tnames)),
#     #                             "logLik", "converged")
#     # # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
#     # attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
#     # attr(out, "REML") <- REML
#     # attr(out, "mincLmerList") <- mincLmerList
#     # attr(out, "atlas") <- attr(anat, "atlas")
#     # attr(out, "definitions") <- attr(anat, "definitions")
#     # attr(out, "anat") <- anat
#     # 
#     # class(out) <- c("anatLmerModel", "anatModel", "matrix")
#     # 
#     # return(out)
#   }
