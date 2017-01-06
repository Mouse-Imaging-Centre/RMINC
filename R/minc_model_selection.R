#' @export
AIC.mincLm <- function(object, ..., k = 2){
  mod_mat <- attr(object, "model")
  dfs <- ncol(mod_mat) + 1 # extra 1 for estimating error scale
  
  k * (dfs - object[,"logLik"])
}

#' @export
AIC.vertexLm <- function(object, ..., k = 2)
  AIC.mincLm(object, ..., k = k)

#' @export
AIC.anatModel <- function(object, ..., k = 2)
  AIC.mincLm(object, ..., k = k)

#' Compute the Corrected AIC for an object
#' 
#' Corrected AIC for finite sample sizes. Generally recommended
#' over AIC as it asymptotically approaches AIC as n->infinity.
#'  
#' @param object An object with an AICc method
#' @param ... additional parameters for methods
#' @details For use with \link{mincLm} objects AICc produces a vector of AICc values, one per voxel.
#' @export
AICc <- function(object, ...)
  UseMethod("AICc")

#' @export
AICc.mincLm <- function(object, ...){
  mod_mat <- attr(object, "model")
  n <- nrow(mod_mat)
  k <- ncol(mod_mat) + 1 # extra 1 for estimating error scale
  
  AIC(object) + 2 * (k + 1) * (k + 2) / (n - k - 2)
}

#' @export
AICc.vertexLm <- function(object, ...)
  AICc.mincLm(object, ...)

#' @export
AICc.anatModel <- function(object, ...)
  AICc.mincLm(object, ...)

#' @export
BIC.mincLm <- function(object, ...){
  mod_mat <- attr(object, "model")
  n <- nrow(mod_mat)
  k <- ncol(mod_mat) + 1 # extra 1 for estimating error scale
  
  -2 * object[,"logLik"] + k * log(n) # omits the - k * log(2*pi) term to match stats::BIC
}

#' @export
BIC.vertexLm <- function(object, ...)
  BIC.mincLm(object, ...)

#' @export
BIC.anatModel <- function(object, ...)
  BIC.mincLm(object, ...)
