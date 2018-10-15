#' Effect Sizes
#'
#' Takes the output of a minc modelling function and computes the unbiased hedges g*
#' and variance of hedges g*
#' @param buffer The results of a vertex/anat/mincLm run
#' @param columns A vector of column names. By default the threshold will
#' be computed for all factor columns.
#' @details This code implements the methods from Nakagawa, S., Cuthill, I.C., 2007. Effect size, confidence interval and statistical significance: a practical guide for biologists. Biol. Rev. Camb. Philos. Soc. 82, 591–605. https://doi.org/10.1111/j.1469-185X.2007.00027.x
#' for computing effect size of group comparisons from a GLM
#' @return A matrix with columns of hedgesg and hedgesg_var for each factor predictor in the GLM
#' or for each column supplied.
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' # read the text file describing the dataset
#' gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
#' # run a linear model relating the data in all voxels to Genotype
#' vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
#' effectsize <- mincEffectSize(vs)
#' }
#' @export
vertexEffectSize <- function(buffer, columns = NULL)
{
  #Nakagawa, S., Cuthill, I.C., 2007. Effect size, confidence interval and statistical significance: a practical guide for biologists. Biol. Rev. Camb. Philos. Soc. 82, 591–605. https://doi.org/10.1111/j.1469-185X.2007.00027.x
  #Unbiased corrector function J from https://en.wikipedia.org/wiki/Effect_size#Hedges'_g
  #Watch out for exploding gamma function
  J <- function(a) {
    out <- tryCatch({
      return (gamma(a / 2) / (sqrt(a / 2) * gamma((a - 1) / 2)))
    },
    error = function(cond) {
      return(1)
    },
    warning = function(cond) {
      return(1)
    })
    return(out)
  }
  # g_biased = t ( n1 + n2 ) / sqrt(n1*n2)*sqrt(df)
  # g_unbiased = g_biased * J(n1 + n2 - 2)
  # g_variance_unbiased = (n1 + n2) / (n1*n2) + g_unbiased**2 / (2 * (n1 + n2 - 2))
  # N = n1+n2
  # NN = n1*n2
  
  originalMincAttrs <- mincAttributes(buffer)
  stattype <- originalMincAttrs$`stat-type`
  # Remove coefficients from buffer
  if (!is.null(stattype)) {
    for (nStat in 1:length(stattype)) {
      if (stattype[nStat] == 'beta' ||
          stattype[nStat] == 'R-squared' ||
          stattype[nStat] == "logLik" || stattype[nStat] == "F") {
        if (!exists('indicesToRemove')) {
          indicesToRemove = nStat
        }
        else {
          indicesToRemove = c(indicesToRemove, nStat)
        }
      }
    }
    if (exists('indicesToRemove')) {
      updatedAttrs <- originalMincAttrs
      updatedAttrs$`stat-type` <-
        updatedAttrs$`stat-type`[-indicesToRemove]
      updatedAttrs$dimnames[[2]] <-
        updatedAttrs$dimnames[[2]][-indicesToRemove]
      
      buffer <- buffer[, -indicesToRemove]
      buffer <- setMincAttributes(buffer, updatedAttrs)
    }
  }
  
  model_call <- attr(buffer, "call")
  model_call[[1]] <- quote(model.matrix)
  names(model_call)[names(model_call) == "formula"] <- ""
  mod_mat <- eval(model_call)
  conts <- attr(mod_mat, "contrasts")
  if (is.null(conts))
    stop("No categorical variables in model, cannot compute g* statistics")
  if (!all(sapply(columns, function(x) any(grepl(x, updatedAttrs$dimnames[[2]])))))
    stop("Columns not found in model")

  if (is.null(columns)) {
    cat_vars <-
      names(conts)[vapply(conts, function(x)
        all(x == "contr.treatment"), logical(1))]
    columns <- lapply(updatedAttrs$dimnames, function(x)
      grep(paste(cat_vars, collapse = "|"), x, value = TRUE))[[2]]
    columns <- gsub("tvalue-", "", columns, fixed = TRUE)
  }
  
  n.cols <- length(columns)
  n.row <- 0
  if (is.matrix(buffer)) {
    n.row <- nrow(buffer)
  }
  else {
    n.row <- length(buffer)
  }
  output <- matrix(1, nrow = n.row, ncol = 2 * n.cols)
  colnames(output) <-
    c(paste0("hedgesg-", columns),
      paste0("hedgesg_var-", columns))
  for (column in columns) {
    n1 <- sum(updatedAttrs$model[, column])
    N <- length(updatedAttrs$model[, column])
    n2 <- N - n1
    NN <- n1 * n2
    df <- updatedAttrs$df[[1]][2]
    output[, paste0("hedgesg-", column)] <-
      buffer[, paste0("tvalue-", column)] * N / (sqrt(NN) * sqrt(df)) * J(N -
                                                                            2)
    output[, paste0("hedgesg_var-", column)] <-
      (N) / (NN) + output[, paste0("hedgesg-", column)] ** 2 / (2 * (N - 2))
  }
  rownames(output) <- rownames(buffer)
  attr(output, "likeVolume") <- attr(buffer, "likeVolume")
  
  # run the garbage collector...
  gcout <- gc()
  
  return(output)
}

#' @describeIn vertexEffectSize mincEffectSize
#' @export
mincEffectSize <- function(buffer, columns = NULL) {
  eff_call <- match.call()
  eff_call[[1]] <- quote(vertexEffectSize)
  eval(eff_call)
}

#' @describeIn vertexEffectSize anatEffectSize
#' @export
anatEffectSize <- function(buffer, columns = NULL) {
  eff_call <- match.call()
  eff_call[[1]] <- quote(vertexEffectSize)
  eval(eff_call)
} 