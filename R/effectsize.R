#' Effect Sizes
#'
#' Takes the output of a minc modelling function and computes the unbiased
#' hedges g* and variance of hedges g*
#' @param buffer The results of a vertex/anat/mincLm run
#' @param predictors A vector of factor predictor names. By default the effect size
#' be computed for all treatment-coded factor columns.
#' @details This code implements the methods from Nakagawa, S., Cuthill, I.C., 2007. Effect size, confidence interval and statistical significance: a practical guide for biologists. Biol. Rev. Camb. Philos. Soc. 82, 591–605. https://doi.org/10.1111/j.1469-185X.2007.00027.x
#' for computing effect size of group comparisons from a GLM.
#' 
#' For now, interactions are explicitly excluded from being predictors. To get
#' effect size for interactions, use the interaction() function to create a new
#' treatment coded factor to use as a predictor.
#' @return A matrix with columns of hedgesg-<factorlevel> and hedgesg_var-<factorlevel> for each factor predictor in the GLM
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
vertexEffectSize <- function(buffer, predictors = NULL)
{
  #Nakagawa, S., Cuthill, I.C., 2007. Effect size, confidence interval and statistical significance: a practical guide for biologists. Biol. Rev. Camb. Philos. Soc. 82, 591–605. https://doi.org/10.1111/j.1469-185X.2007.00027.x
  #Original unbiased corrector from paper replaced with
  #unbiased corrector function J from https://en.wikipedia.org/wiki/Effect_size#Hedges'_g
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
  # n1 study group
  # n2 reference group
  # N = n1+n2
  # NN = n1*n2

  #Scrub columns not needed
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

  #Extract model, contrasts from model call
  model_call <- attr(buffer, "call")
  model_call[[1]] <- quote(model.matrix)
  names(model_call)[names(model_call) == "formula"] <- ""
  mod_mat <- eval(model_call)
  conts <- attr(mod_mat, "contrasts")
  cat_vars <-
    names(conts)[vapply(conts, function(x)
      all(x == "contr.treatment"), logical(1))]
  
  #Checking for assumptions regarding computing
  if (is.null(conts))
    stop("No factors in model, cannot compute g* statistics")
  
  if (is.null(cat_vars))
    stop("No treatment-coded factors in model, cannot compute g* statistics")
  
  if (!is.null(predictors) &&
      !all(predictors %in% colnames(updatedAttrs$model)))
    stop("Supplied predictors not found in model")
  
  if (!is.null(predictors) &&
      !all(grepl(paste(cat_vars, collapse = "|"), predictors)))
    stop("Supplied predictors are not treatment contrasts")
  
  if (!is.null(predictors) && any(grepl(":", predictors)))
      stop("Interactions in predictors are not currently supported, 
           generate treatment contrasts using interaction()")

  #If predictors aren't provided, ennumerate factors from the model
  if (is.null(predictors)) {
    predictors <- grep(":", grep(paste(cat_vars, collapse = "|"), updatedAttrs$dimnames[[2]], value = TRUE), invert = TRUE, value = TRUE)
    predictors <- gsub("tvalue-", "", predictors, fixed = TRUE)
  }

  
  #Compute reference group size, find number of total subjects, subtract all
  #instaces of the treatment-coded factor, remaining number is reference group
  referencegroupsize <-
    matrix(1, nrow = 1, ncol = length(predictors))
  colnames(referencegroupsize) <- predictors

  for (var in cat_vars) {
    referencegroupsize[,  grep(var, predictors, value = TRUE)] <-
      NROW(updatedAttrs$model[, grep(":", grep(var, colnames(updatedAttrs$model), value = TRUE), invert = TRUE, value = TRUE)]) -
      sum(updatedAttrs$model[, grep(":", grep(var, colnames(updatedAttrs$model), value = TRUE), invert = TRUE, value = TRUE)])
  }

  n.cols <- length(predictors)
  n.row <- 0
  if (is.matrix(buffer)) {
    n.row <- nrow(buffer)
  }
  else {
    n.row <- length(buffer)
  }
  output <- matrix(1, nrow = n.row, ncol = 2 * n.cols)
  colnames(output) <-
    c(paste0("hedgesg-", predictors),
      paste0("hedgesg_var-", predictors))
  for (column in predictors) {
    n1 <- sum(updatedAttrs$model[, column])
    n2 <- referencegroupsize[, column]
    N <- n1 + n2
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
mincEffectSize <- function(buffer, predictors = NULL) {
  eff_call <- match.call()
  eff_call[[1]] <- quote(vertexEffectSize)
  eval(eff_call)
}

#' @describeIn vertexEffectSize anatEffectSize
#' @export
anatEffectSize <- function(buffer, predictors = NULL) {
  eff_call <- match.call()
  eff_call[[1]] <- quote(vertexEffectSize)
  eval(eff_call)
}
