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
vertexMean <- function(filenames) 
{
  vertexData = vertexTable(filenames)
  return(rowMeans(vertexData))
  
} 

#' @describeIn vertexSummaries sum
#' @export
vertexSum <- function(filenames) 
{
  vertexData = vertexTable(filenames)
  return(rowSums(vertexData))
  
} 

#' @describeIn vertexSummaries var
#' @export
vertexVar <- function(filenames) 
{
  vertexData = vertexTable(filenames)
  return(apply(vertexData,1,var))
  
} 

#' @describeIn vertexSummaries standard deviation
#' @export
vertexSd<- function(filenames) 
{
  vertexData = vertexTable(filenames)
  return(apply(vertexData,1,sd))
  
} 

#' @description This function is used to compute an arbitrary function of every region in a set of vertex files.
#' @name vertexApply
#' @title Apply function over vertex Files
#' @param filenames vertex file names
#' @param function.string The function which to apply. Can only take a single
#' argument, which has to be 'x'.
#' @return  out: The output will be a single vector containing as many
#'          elements as there are vertices in the input variable. 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- vertexApply(gf$CIVETFILES$nativeRMStlink20mmleft,quote(mean(x)))
#' }
#' @export
vertexApply <- function(filenames,function.string) 
{
  # Load the data
  vertexData = vertexTable(filenames)
  
  # In order to maintain the same interface as mincApply, the (x) part needs to be stripped
  function.string = gsub('(x)','', function.string, fixed = TRUE)
  
  # The apply part (transpose to match output of mincApply)
  results <- t(apply(vertexData,1,function.string[1]))
  
  attr(results, "likeFile") <- filenames[1]
  
  # run the garbage collector...
  gcout <- gc()
  
  return(results)
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
vertexAnova <- function(formula, data, subset=NULL) {
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
  data.matrix <- vertexTable(filenames)
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
vertexLm <- function(formula, data, subset=NULL) {
  
  # Build model.frame
  m <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  
  
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
    parseLmOutput$data.matrix.left  <- vertexTable(parseLmOutput$data.matrix.left)
    parseLmOutput$data.matrix.right <- vertexTable(parseLmOutput$data.matrix.right)
  }
  else {
    filenames <- as.character(mf[,1])
    parseLmOutput$mmatrix <- model.matrix(formula, mf)	
    parseLmOutput$data.matrix.left <- vertexTable(filenames)
    parseLmOutput$rows = colnames(parseLmOutput$mmatrix)
  } 
  
  result <- .Call("vertex_lm_loop",
                  parseLmOutput$data.matrix.left,
                  parseLmOutput$data.matrix.right,
                  parseLmOutput$mmatrix,
                  PACKAGE="RMINC") 
  
  attr(result, "likeVolume") <- as.character(mf[,1])[1]
  attr(result, "model") <- as.matrix(parseLmOutput$mmatrix)
  attr(result, "filenames") <- as.character(mf[,1])
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))
  
  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2+1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  betaNames = paste('beta-', parseLmOutput$rows, sep='')
  tnames = paste('tvalue-', parseLmOutput$rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("vertexMultiDim", "matrix")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

vertexApplyRCPP <-
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
  

vertexLmer <-
  function(formula, data, mask=NULL, parallel=NULL,
           REML=TRUE, control=lmerControl(), start=NULL, 
           verbose=0L, temp_dir = getwd(), safely = FALSE, 
           cleanup = TRUE) {
    
  }