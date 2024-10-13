#' flexLm
#'
#' @param formula 
#' @param data 
#' @param y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
flexLm <- function(formula, data, y, ...) {
  # collect all the datavars that were given in the ... flexible argument list
  dataVars <- list(...)
  
  # they have to be named; throw an error if not
  if(any(names(dataVars)=="") | any(is.null(names(dataVars)))) {
    stop("All data variables have to be named")
  }
  # could improve this by varnames=lapply(substitute(list(...))[-1], deparse) to get variable names used
  
  # dimensions of y and all dataVars have to be the same
  if (length(dataVars) > 0) {
    for (i in 1:length(dataVars)) {
      if (any(dim(dataVars[[i]]) != dim(y))) {
        stop("All data variables and the y variable have to be the same dimension")
      }
    }
  }  
  
  # now we set up the model matrix, using fake data in the first instance
  fakeData <- data
  fakeData$y <- y[1,]
  fakeData2 <- fakeData
  if (length(dataVars) > 0) {
    for (i in 1:length(dataVars)) {
      # two sets of fake data - will let us evaluate whether the formula shortcut is working
      fakeData[,names(dataVars)[i]] <- rnorm(nrow(fakeData))
      fakeData2[,names(dataVars)[i]] <- 0
    }
  }
  
  # generate two model matrices, one from each set of fake data
  mmatrix <- model.matrix(formula(formula), fakeData)
  mmatrix2 <- model.matrix(formula(formula), fakeData2)
  
  # find the matrix columns that differ between the model matrices - these
  # will be the ones that depend on data that varies at every voxel/vertex/whatever
  mmatrixDiff <- mmatrix
  mmatrixDiff[,] <- as.integer( !(mmatrix == mmatrix2))
  
  # now test wheter a simple multiplication of the two matrices produces the correct model matrix
  # if it does we can be much faster in looping over voxels/vertices/whatever, whereas
  # if it doesn't we have to evaluate the formula every time.
  # the shortcut should work if there are no transformations of the voxel data (i.e. taking squares, splines, etc)
  
  shortcutWorks <- all( mmatrix == (mmatrix2 + (mmatrixDiff * as.vector(fakeData[,names(dataVars)][[1]]) )))
  
  # create the output matrix
  out <- matrix(nrow=nrow(y), ncol=ncol(mmatrix)*2)
  
  # now loop over every voxel/vertex/whatever
  for (i in 1:nrow(y)) {
    
    # if the matrix multiplication shortcut works use it
    if (shortcutWorks) {
      tmp <- fastLm(mmatrix2 + (mmatrixDiff * dataVars[[1]][i,]), y[i,])
    } else {
      # otherwise re-evaluate the formula at every voxel after assign the right
      # values to the data frame
      if (length(dataVars) > 0) {
        for (j in 1:length(dataVars)) {
          data[,names(dataVars)[j]] <- dataVars[[j]][i,]
        }
      }
      mmatrix <- model.matrix(formula(formula), data)
      tmp <- fastLm(mmatrix, y[i,])
    }
    # keep the coefficients and t statistics
    out[i,] <- c(tmp$coefficients, tmp$coefficients/tmp$stderr)
  }
  
  # assign names to the outputs and return
  nnames <- c(names(tmp$coefficients), paste("tstatistic", names(tmp$coefficients),  sep=":"))
  colnames(out) <- nnames
  return(out)
}