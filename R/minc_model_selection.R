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


check_model_list <- function(model_list, model_type){
  ## Function to check model_list represents a valid set of minc models
  
  if(length(model_list) < 2)
    stop("more than one model must be supplied")
  
  if(any(!sapply(model_list, inherits, model_type)))
    stop("all models must be the same type")
  
  model_lengths <- sapply(model_list, nrow)
  
  if(any(model_lengths != model_lengths[1]))
    stop("all models must have the same number of rows")
}

#' Compare a set of massively-univariate models
#' 
#' For each response (voxel, vertex, or structure) compare
#' a set of linear model formulations with the criterion of your
#' choice (e.g. AIC, AICc, BIC). 
#' 
#' @param object A model or a list of models, typically 
#' \link{mincLm}, \link{vertexLm}, or \link{anatLm} results.
#' @param ... additional models
#' @param metric A function to apply to the models that extracts a result
#' for each independent sub-model. Typical choices are \link{AIC}, \link{AICc}, 
#' and \link{BIC}. Please note that metrics are considered such that lower is 
#' better (in following AIC). To use a positive metric create a wrapper function
#' that performs the negation, for example, to use the un-modified log-likelihood 
#' you could pass \code{ metric = function(minc_model){ -minc_model[,"logLik"]} }
#' @return A sub-model x n models \code{model_comparison} matrix with the metric of interest.
#' @export
compare_models <- function(object, ..., metric = AICc)
  UseMethod("compare_models")

#' @export
compare_models.mincLm <- function(object, ..., metric = AICc){
  models <- list(object, ...)
  metric <- match.fun(metric)
  
  check_model_list(models, class(object))
  
  metric_table <-
    lapply(models, metric) %>%
    Reduce(cbind, ., NULL)
  
  attr(metric_table, "class")    <- c("model_comparison", "matrix")
  attr(metric_table, "formulae") <- lapply(models, function(m) attr(m, "call")[["formula"]])
  attr(metric_table, "metric")   <- deparse(substitute(metric))
  
  metric_table
}

#' @export
compare_models.vertexLm <- function(object, ..., metric = AICc)
  compare_models.mincLm(object, ..., metric = metric)

#' @export
compare_models.anatModel <- function(object, ..., metric = AICc)
  compare_models.mincLm(object, ..., metric = metric)

#' @export
print.model_comparison <- function(x, n = min(6, nrow(x)), width = min(6, ncol(x)), ...){
  formula_strings <- sapply(attr(x, "formulae")[seq_len(width)], deparse)
  
  cat("model_comparison matrix, showing ", n, " of ", nrow(x), " rows, and ", width, " of ", ncol(x), "columns\nprint more by supplying n and width to the print function\n")
  cat("showing models:\n")
  cat(paste0(seq_along(formula_strings), ": ", formula_strings), sep = "\n")
  cat("\n")
  
  sub_x <- x[seq_len(n), seq_len(width)]
  
  print(sub_x)
  invisible(x)
}

#' @export
summary.model_comparison <- function(object, ...){
  victors <- unlist(apply(object, 1, function(x) which(x == min(x))))
  wins_table <- 
    table(victors) %>% 
    { data_frame(model  = names(.), wins = .)}
  
  formula_strings <- sapply(attr(object, "formulae"), deparse)
  formula_table <- 
    data_frame(formula = formula_strings, model = as.character(seq_along(formula_strings)))
  
  wins_frame <-
    full_join(formula_table, wins_table, by = "model") %>%
    setNA(0) %>%
    select(.data$model, .data$formula, .data$wins) %>%
    arrange(.data$wins)
  
  wins_frame
}
