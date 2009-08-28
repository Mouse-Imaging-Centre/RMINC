###            Fit a general linear mixed effects model

library('nlme')

newlme <-
  ## fits general linear mixed effects model by maximum likelihood, or
  ## residual maximum likelihood using Newton-Raphson algorithm.
  function(fixed,
	   data = sys.frame(sys.parent()),
	   random,
	   correlation = NULL,
	   weights = NULL,
	   subset,
	   method = c("REML", "ML"),
	   na.action = na.fail,
	   control = list(),
           contrasts = NULL, keep.data = TRUE)
  UseMethod("newlme")

newlme.formula <-
  function(fixed,
	   data = sys.frame(sys.parent()),
	   random = pdSymm( eval( as.call( fixed[ -2 ] ) ) ),
	   correlation = NULL,
	   weights = NULL,
	   subset,
	   method = c("REML", "ML"),
	   na.action = na.fail,
	   control = list(),
           contrasts = NULL,
           keep.data = TRUE)
{   
  Call <- match.call()
  miss.data <- missing(data) || !is.data.frame(data)

  ## control parameters
  controlvals <- newlmeControl()

  ##
  ## checking arguments
  ##
  if (!inherits(fixed, "formula") || length(fixed) != 3) {
    stop("\nFixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  method <- match.arg(method)
  REML <- method == "REML"
  reSt <- reStruct(random, REML = REML, data = NULL)
  groups <- newgetGroupsFormula(reSt)

  ## create an lme structure containing the random effects model and plug-ins

  lmeSt <- newlmeStruct(reStruct = reSt, corStruct = correlation,
		     varStruct = newvarFunc(weights))

  ## extract a data frame with enough information to evaluate
  ## fixed, groups, reStruct, corStruct, and varStruct
  mfArgs <- list(formula = newasOneFormula(formula(lmeSt), fixed, groups),
		 data = data, na.action = na.action)
  
  mfArgs$drop.unused.levels <- TRUE
  dataMix <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMix)	# preserve the original order
  for(i in names(contrasts))            # handle contrasts statement
      contrasts(dataMix[[i]]) = contrasts[[i]]

  ## sort the model.frame by groups and get the matrices and parameters
  ## used in the estimation procedures

  grps <- newgetGroups(dataMix, groups)  

  ## ordering data by groups
  if (inherits(grps, "factor")) {	# single level
    ord <- order(grps)	#"order" treats a single named argument peculiarly
    grps <- data.frame(grps)
    row.names(grps) <- origOrder
    names(grps) <- as.character(deparse((groups[[2]])))
  } else {
    ord <- do.call("order", grps)
    ## making group levels unique
    for(i in 2:ncol(grps)) {
      grps[, i] <-
        as.factor(paste(as.character(grps[, i-1]), as.character(grps[,i]),
                        sep = "/"))
      NULL
    }
  }
  
  grps <- grps[ord, , drop = FALSE]
  dataMix <- dataMix[ord, ,drop = FALSE]
  #revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
  ## obtaining basic model matrices
  N <- nrow(grps)
  Z <- model.matrix(reSt, dataMix)
  ncols <- attr(Z, "ncols")
  Names(lmeSt$reStruct) <- attr(Z, "nams")
  ## keeping the contrasts for later use in predict

  contr <- attr(Z, "contr")
  X <- model.frame(fixed, dataMix)
  #Terms <- attr(X, "terms")
  auxContr <- lapply(X, function(el)
		     if (inherits(el, "factor") &&
                         length(levels(el)) > 1) contrasts(el))
  contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
  contr <- contr[!unlist(lapply(contr, is.null))]
  X <- model.matrix(fixed, data=X)

  ncols <- c(ncols, dim(X)[2], 1)
  Q <- ncol(grps)

  dims <- newMEdims(grps, ncols)

  listNncols <- c(N, sum(ncols))
  listrownames <- list(row.names(dataMix), c(colnames(Z), colnames(X),deparse(fixed[[2]])))

  attr(lmeSt, "conLin") <-
    list(Xy = array(c(Z, X, eval(fixed[[2]], dataMix)), listNncols, listrownames),
	 dims = dims, logLik = 0)
  lmeSt <- newInitialize(lmeSt, dataMix, grps, control = controlvals)

  dataMix <<-dataMix
  X<<-X
  Z<<-Z
  grps<<-grps
  lmeSt<<-lmeSt
  controlvals<<-controlvals
  dims <<- dims
  listNncols <<- listNncols
  listrownames <<- listrownames

  storedTstats <- newlmeLoop(dataMix, X, Z, grps, lmeSt, controlvals, dims, listNncols, listrownames, eval(fixed[[2]], dataMix))
  return(storedTstats)
}

newlmeLoop <- function(dataMix, X, Z, grps, lmeSt, controlvals, dims, listNncols, listrownames, y)
{
  #y <- eval(fixed[[2]], dataMix)
  #dataMix$x <-  y[attr(dataMix, "row.names")]

  y <-  y[attr(dataMix, "row.names")]

  ## checking if enough observations per group to estimate ranef
  #tmpDims <- attr(lmeSt, "conLin")$dims

  ## degrees of freedom for testing fixed effects
#   fixDF <- newgetFixDF(X, grps, attr(lmeSt, "conLin")$dims$ngrps,
#                     terms = Terms) 
  
  ## initialization
#   lmeSt <- newInitialize(lmeSt, dataMix, grps, control = controlvals)
#   parMap <- attr(lmeSt, "pmap")

  ## creating the condensed linear model  
  attr(lmeSt, "conLin") <-
    list(Xy = array(c(Z, X, y), listNncols, listrownames),
	 dims = dims, logLik = 0)
  #tmpDims <- attr(lmeSt, "conLin")$dims

  ## Checking possibility of single decomposition
#   if (length(lmeSt) == 1)  {	# reStruct only, can do one decomposition
#     ## need to save conLin for calculating fitted values and residuals
# #     oldConLin <- attr(lmeSt, "conLin")
#     decomp <- TRUE
#     attr(lmeSt, "conLin") <- newMEdecomp(attr(lmeSt, "conLin"))
#   } else decomp <- FALSE
  ##
  ## getting the linear mixed effects fit object,
  ## possibly iterating for variance functions
  ##

  numIter <- 0
  repeat {
    oldPars <- coef(lmeSt)
    optRes <- if (controlvals$opt == "nlminb") {
        nlminb(c(coef(lmeSt)),
               function(lmePars) -logLik(lmeSt, lmePars),
               control = list(iter.max = controlvals$msMaxIter,
               eval.max = controlvals$msMaxEval,
               trace = controlvals$msVerbose))
    } else {
        optim(c(coef(lmeSt)),
              function(lmePars) -logLik(lmeSt, lmePars),
              control = list(trace = controlvals$msVerbose,
              maxit = controlvals$msMaxIter,
              reltol = if(numIter == 0) controlvals$msTol
              else 100*.Machine$double.eps),
              method = controlvals$optimMethod)
    }
    numIter0 <- NULL
    coef(lmeSt) <- optRes$par
    attr(lmeSt, "lmeFit") <- newMEestimate(lmeSt, grps)
    ## checking if any updating is needed
    if (!needUpdate(lmeSt)) {
	if (optRes$convergence) {
	    msg <- paste(controlvals$opt, " problem, convergence error code = ",
			 optRes$convergence, "\n  message = ", optRes$message,
			 sep='')
	    if(!controlvals$returnObject) {
		#stop(msg)
                warning(msg)
            }
	    else
		warning(msg)
	}
	break
    }

    ## updating the fit information
    numIter <- numIter + 1
    lmeSt <- update(lmeSt, dataMix)
    ## calculating the convergence criterion
    aConv <- coef(lmeSt)
    conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
    aConv <- NULL
    for(i in names(lmeSt)) {
	if (any(parMap[,i])) {
	    aConv <- c(aConv, max(conv[parMap[,i]]))
	    names(aConv)[length(aConv)] <- i
	}
    }
    if (max(aConv) <= controlvals$tolerance) {
	break
    }
    if (numIter > controlvals$maxIter) {
	msg <- paste("Maximum number of iterations",
		     "(lmeControl(maxIter)) reached without convergence.")
	if (controlvals$returnObject) {
	    warning(msg)
	    break
	} else
	    stop(msg)
    }

  } ## end{repeat}

  ## wrapping up
  lmeFit <- attr(lmeSt, "lmeFit")

  #attr(fixDF, "varFixFact") <- varFix <- lmeFit$sigma * lmeFit$varFix
  varFix <- crossprod(lmeFit$sigma * lmeFit$varFix)

  names(lmeFit$beta) <- colnames(X)
  dimnames(varFix) <- list(colnames(X), colnames(X))
  
  ##
  ## fitted.values and residuals (in original order)
  ##
#   Fitted <- fitted(lmeSt, level = 0:Q,
# 		   conLin = if (decomp) oldConLin else attr(lmeSt, "conLin"))[
# 		   revOrder, , drop = FALSE]
#   Resid <- y[revOrder] - Fitted
#   rownames(Resid) <- rownames(Fitted) <- origOrder
#   attr(Resid, "std") <- lmeFit$sigma/(newvarWeights(lmeSt)[revOrder])
  ## putting groups back in original order
#   grps <- grps[revOrder, , drop = FALSE]
  ## making random effects estimates consistently ordered
#  for(i in names(lmeSt$reStruct)) {
#    lmeFit$b[[i]] <- lmeFit$b[[i]][unique(as.character(grps[, i])),, drop = F]
#    NULL
#  }
  ## inverting back reStruct
#   lmeSt$reStruct <- solve(lmeSt$reStruct)
  ## saving part of dims
#   dims <- attr(lmeSt, "conLin")$dims[c("N", "Q", "qvec", "ngrps", "ncol")]
  ## getting the approximate var-cov of the parameters
#   if (controlvals$apVar) {
#     apVar <- newlmeApVar(lmeSt, lmeFit$sigma,
# 		      .relStep = controlvals[[".relStep"]],
#                       minAbsPar = controlvals[["minAbsParApVar"]],
# 		      natural = controlvals[["natural"]])
#   } else {
#     apVar <- "Approximate variance-covariance matrix not available"
#   }
#   ## getting rid of condensed linear model and fit
#   attr(lmeSt, "conLin") <- NULL
#   attr(lmeSt, "lmeFit") <- NULL
#   ##
#   ## creating the lme object
#   ##
#   estOut <- list(modelStruct = lmeSt,
# 		 dims = dims,
# 		 contrasts = contr,
# 		 coefficients = list(
# 		     fixed = lmeFit$beta,
# 		     random = lmeFit$b),
# 		 varFix = varFix,
# 		 sigma = lmeFit$sigma,
# 		 apVar = 1,#apVar,
# 		 logLik = lmeFit$logLik,
# 		 numIter = if (needUpdate(lmeSt)) numIter
# 		   else numIter0,
# 		 groups = grps,
# 		 call = Call,
#                  terms = Terms,
# 		 method = method,
# 		 fitted = Fitted,
# 		 residuals = Resid,
#                  fixDF = fixDF,
#                  na.action = attr(dataMix, "na.action"))
# #   if (keep.data && !miss.data) estOut$data <- data
# #   if (inherits(data, "groupedData")) {
# #     ## saving labels and units for plots
# #     attr(estOut, "units") <- attr(data, "units")
# #     attr(estOut, "labels") <- attr(data, "labels")
# #   }
#   class(estOut) <- "lme"
#   print(estOut)

  storedTstats <- matrix(nrow=1, ncol= length(lmeFit$beta))
  colnames(storedTstats) <- names(lmeFit$beta)

  storedTstats[1,] <- as.matrix(lmeFit$beta)/sqrt(diag(as.matrix(varFix)))

  return(storedTstats)
}

### Auxiliary functions used internally in lme and its methods

newgetFixDF <-
  function(X, grps, ngrps, assign = attr(X, "assign"), terms)
{
  ## calculates degrees of freedom for fixed effects Wald tests
  if (!is.list(assign)) {               # in R
    namTerms <- attr(terms, "term.labels")
    if (attr(terms, "intercept") > 0) {
      namTerms <- c("(Intercept)", namTerms)
    }
    namTerms <- factor(assign, labels = namTerms)
    assign <- split(order(assign), namTerms)
  }
  ## function to check if a vector is (nearly) a multiple of (1,1,...,1)
  const <- function(x, tolerance = sqrt(.Machine$double.eps)) {
      if (length(x) < 1) return(NA)
      x <- as.numeric(x)
      if (x[1] == 0.) return(all(abs(x) < tolerance))
      all(abs((x/x[1] - 1.)) < tolerance)
  }
  N <- nrow(X)
  p <- ncol(X)
  Q <- ncol(grps)
  Qp1 <- Q + 1
  namX <- colnames(X)
  ngrps <- rev(ngrps)[-(1:2)]
  stratNam <- c(names(ngrps), "Residual")
  dfX <- dfTerms <- c(ngrps, N) - c(0, ngrps)
  names(dfX) <- names(dfTerms) <- stratNam
  valX <- double(p)
  names(valX) <- namX
  namTerms <- names(assign)
  valTerms <- double(length(assign))
  names(valTerms) <- namTerms
  if (any(notIntX <- !apply(X, 2, const))) {
      ## percentage of groups for which columns of X are inner
      innP <- array(c(rep(1, p),
                      .C("inner_perc_table",
                         as.double(X),
                         as.integer(unlist(grps)),
                         as.integer(p),
                         as.integer(Q),
                         as.integer(N),
                         val = double(p * Q))[["val"]]), c(p, Qp1),
                    list(namX, stratNam))
    ## strata in which columns of X are estimated
    ## ignoring fractional inner percentages for now
    stratX <- stratNam[apply(innP, 1, function(el, index) max(index[el > 0]),
                             index = 1:Qp1)]
    ## strata in which terms are estimated
    notIntTerms <- unlist(lapply(assign,
                                 function(el, notIntX) {
                                   any(notIntX[el])
                                 }, notIntX = notIntX))
    stratTerms <- stratNam[unlist(lapply(assign,
                          function(el, stratX, stratNam) {
                            max(match(stratX[el], stratNam))
                          },
                       stratX = stratX, stratNam = stratNam))][notIntTerms]
    stratX <- stratX[notIntX]
    xDF <- table(stratX)
    dfX[names(xDF)] <- dfX[names(xDF)] - xDF
    if (!all(notIntX)) {                # correcting df for intercept
      dfX[1] <- dfX[1] - 1
    } else {
      dfX[-1] <- dfX[-1] + 1
    }
    valX[notIntX] <- dfX[stratX]
    ## number of parameters in each term
    pTerms <- unlist(lapply(assign, length))[notIntTerms]
    tDF <- tapply(pTerms, stratTerms, sum)
    dfTerms[names(tDF)] <- dfTerms[names(tDF)] - tDF
    if (!all(notIntTerms)) {
      dfTerms[1] <- dfTerms[1] - 1
    } else {
      dfTerms[-1] <- dfTerms[-1] + 1
    }
    valTerms[notIntTerms] <- dfTerms[stratTerms]
  } else {
    notIntTerms <- unlist(lapply(assign,
                                 function(el, notIntX) {
                                   any(notIntX[el])
                                 }, notIntX = notIntX))
  }
  if (!all(notIntX)) {  #intercept included
    valX[!notIntX] <- max(dfX)
    if (!all(notIntTerms)) {
      valTerms[!notIntTerms] <- max(dfTerms)
    }
  }
  val <- list(X = valX, terms = valTerms)
  attr(val, "assign") <- assign
  val
}

newlmeApVar <-
  function(lmeSt, sigma, conLin = attr(lmeSt, "conLin"),
           .relStep = (.Machine$double.eps)^(1/3), minAbsPar = 0,
           natural = TRUE)
{
  ## calculate approximate variance-covariance matrix of all parameters
  ## except the fixed effects. By default, uses natural parametrization for
  ## for pdSymm matrices
  fullLmeLogLik <-
    function(Pars, object, conLin, dims, N, settings) {
      ## logLik as a function of sigma and coef(lmeSt)
      npar <- length(Pars)
      sigma <- exp(Pars[npar])              # within-group std. dev.
      Pars <- Pars[-npar]
      coef(object) <- Pars
      if ((lO <- length(object)) > 1) {
	for(i in lO:2) {
	  conLin <- newrecalc(object[[i]], conLin)
	  NULL
	}
      }
      val <- .C("mixed_loglik",
		as.double(conLin$Xy),
		as.integer(unlist(dims)),
		as.double(sigma * unlist(pdFactor(solve(object$reStruct)))),
		as.integer(settings),
		logLik = double(1),
		lRSS = double(1))[c("logLik", "lRSS")]
      aux <- (exp(val[["lRSS"]])/sigma)^2
      conLin[["logLik"]] + val[["logLik"]] + (N * log(aux) - aux)/2
    }
  dims <- conLin$dims
  sett <- attr(lmeSt, "settings")
  N <- dims$N - sett[1] * dims$ncol[dims$Q + 1]
  sett[2:3] <- c(1, 0)			# asDelta = TRUE and no grad/Hess
  conLin[["logLik"]] <- 0               # making sure
  sig2 <- sigma * sigma
  reSt <- lmeSt[["reStruct"]]
  for(i in seq_along(reSt)) {
    matrix(reSt[[i]]) <- as.double(sig2) * pdMatrix(reSt[[i]])
    if (inherits(reSt[[i]], "pdSymm") && natural) {
      reSt[[i]] <- pdNatural(reSt[[i]])
    }
    if (inherits(reSt[[i]], "pdBlocked") && natural) {
      for(j in seq_along(reSt[[i]])) {
        if (inherits(reSt[[i]][[j]], "pdSymm")) {
          reSt[[i]][[j]] <- pdNatural(reSt[[i]][[j]])
        }
      }
    }
  }
  lmeSt[["reStruct"]] <- reSt
  cSt <- lmeSt[["corStruct"]]
  if (!is.null(cSt) && inherits(cSt, "corSymm") && natural) {
    cStNatPar <- coef(cSt, unconstrained = FALSE)
    class(cSt) <- c("corNatural", "corStruct")
    coef(cSt) <- log((cStNatPar + 1)/(1 - cStNatPar))
    lmeSt[["corStruct"]] <- cSt
  }
  Pars <- c(coef(lmeSt), lSigma = log(sigma))
  val <- newfdHess(Pars, fullLmeLogLik, lmeSt, conLin, dims, N, sett,
		.relStep = .relStep, minAbsPar = minAbsPar)[["Hessian"]]
  if (all(eigen(val)$values < 0)) {
    ## negative definite - OK
    val <- solve(-val)
    nP <- names(Pars)
    dimnames(val) <- list(nP, nP)
    attr(val, "Pars") <- Pars
    attr(val, "natural") <- natural
    val
  } else {
    ## problem - solution is not a maximum
    "Non-positive definite approximate variance-covariance"
  }
}

newMEdecomp <-
 function(conLin)
  ## decompose a condensed linear model.  Returns another condensed
  ## linear model
{
  dims <- conLin$dims
  if (dims[["StrRows"]] >= dims[["ZXrows"]]) {
    ## no pint in doing the decomposition
    return(conLin)
  }
  dc <- array(.C("mixed_decomp",
		 as.double(conLin$Xy),
		 as.integer(unlist(dims)))[[1]],
	      c(dims$StrRows, dims$ZXcols))
  dims$ZXrows <- dims$StrRows
  dims$ZXoff <- dims$DecOff
  dims$ZXlen <- dims$DecLen
  conLin[c("Xy", "dims")] <- list(Xy = dc, dims = dims)
  conLin
}

newMEEM <-
  function(object, conLin, niter = 0)
  ## perform niter iterations of the EM algorithm for conLin
  ## assumes that object is in precision form
{
  if (niter > 0) {
    dd <- conLin$dims
    pdCl <- attr(object, "settings")[-(1:3)]
    pdCl[pdCl == -1] <- 0
    precvec <- unlist(pdFactor(object))
    zz <- .C("mixed_EM",
	     as.double(conLin$Xy),
	     as.integer(unlist(dd)),
	     precvec = as.double(precvec),
	     as.integer(niter),
	     as.integer(pdCl),
	     as.integer(attr(object, "settings")[1]),
	     double(1),
	     double(length(precvec)),
	     double(1))[["precvec"]]
    Prec <- vector("list", length(object))
    names(Prec) <- names(object)
    for (i in seq_along(object)) {
      len <- dd$qvec[i]^2
      matrix(object[[i]]) <-
        crossprod(matrix(zz[1:len + dd$DmOff[i]], ncol = dd$qvec[i]))
    }
  }
  object
}

newMEestimate <-
  function(object, groups, conLin = attr(object, "conLin"))
{
  dd <- conLin$dims
  nc <- dd$ncol
  REML <- attr(object$reStruct, "settings")[1]
  Q <- dd$Q
  rConLin <- newrecalc(object, conLin)
  zz <- .C("mixed_estimate",
	   as.double(rConLin$Xy),
	   as.integer(unlist(dd)),
	   as.double(unlist(pdFactor(object$reStruct))),
	   as.integer(REML),
	   double(1),
	   estimates = double(dd$StrRows * dd$ZXcols),
	   as.logical(FALSE))[["estimates"]]
  estimates <- array(zz, c(dd$StrRows, dd$ZXcols))
  resp <- estimates[ , dd$ZXcols]
  reSt <- object$reStruct
  nam <- names(reSt)
  val <- vector(mode = "list", length = Q)
  names(val) <- nam
  start <- dd$StrRows * c(0, cumsum(nc))
  for (i in seq_along(reSt)) {
    val[[i]] <-
      matrix(resp[as.vector(outer(1:(nc[i]), dd$SToff[[i]] - start[i], "+"))],
	     ncol = nc[i], byrow = TRUE,
	     dimnames = list(unique(as.character(groups[, nam[i]])),
		 Names(reSt[[i]])))
    NULL
  }
  p <- nc[Q + 1]
  N <- dd$N - REML * p
  dimE <- dim(estimates)
  list(logLik = N * (log(N) - (1 + log(2 * pi)))/2 + rConLin$logLik,
       b = rev(val),
       beta = resp[dimE[1] - (p:1)],
       sigma = abs(resp[dimE[1]])/sqrt(N),
       varFix = t(solve(estimates[dimE[1]-(p:1), dimE[2]-(p:1), drop = FALSE])))
}

newMEdims <-
  function(groups, ncols)
{
  ## define constants used in matrix decompositions and log-lik calculations
  ## first need some local functions
  lengths <-
    ## returns the group lengths from a vector of last rows in the group
    function(lstrow) diff(c(0, lstrow))
  offsets <-
    ## converts total number of columns(N), columns per level(ncols), and
    ## a list of group lengths to offsets in C arrays
    function(N, ncols, lstrow, triangle = FALSE)
  {
    pop <- function(x) x[-length(x)]
    cstart <- c(0, cumsum(N * ncols))
    for (i in seq_along(lstrow)) {
      lstrow[[i]] <- cstart[i] +
        if (triangle) {
          lstrow[[i]] - ncols[i]        # storage offsets style
        } else {
          pop(c(0, lstrow[[i]]))        # decomposition style
        }
    }
    lstrow
  }
  Q <- ncol(groups)                     # number of levels
  N <- nrow(groups)                     # number of observations
  ## 'isLast' indicates if the row is the last row in the group at that level.
  ## this version propagates changes from outer groups to inner groups
#  isLast <- (array(unlist(lapply(c(rev(as.list(groups)),
#                                 list(X = rep(0, N), y = rep(0, N))),
#                                function(x) c(0 != diff(codes(x)), TRUE))),
#                  c(N, Q+2), list(NULL, c(rev(names(groups)), "X", "y")))
#             %*% (row(diag(Q+2)) >= col(diag(Q+2)))) != 0
  ## this version does not propagate changes from outer to inner.
  isLast <- array(FALSE, dim(groups) + c(0, 2),
                  list(NULL, c(rev(names(groups)), "X", "y")))
  for(i in 1:Q) {
    isLast[, Q + 1 - i] <- c(0 != diff(as.integer(groups[[i]])), TRUE)
  }
  isLast[N,  ] <- TRUE
  lastRow <- apply(isLast, 2, function(x) seq_along(x)[x])
  if(!is.list(lastRow)) {
    nm <- names(lastRow)
    lastRow <- as.list(lastRow)
    names(lastRow) <- nm
  }

  isLast <- t(isLast)
  strSizes <- cumsum(ncols * isLast) * isLast # required storage sizes
  lastStr <- apply(t(strSizes), 2, function(x) x[x != 0])
  if(!is.list(lastStr)) {
    nm <- names(lastStr)
    lastStr <- as.list(lastStr)
    names(lastStr) <- nm
  }
  strRows <- max(lastStr[[length(lastStr)]])
  lastBlock <- vector("list", Q)
  names(lastBlock) <- rownames(strSizes)[1:Q]
  for(i in 1:Q) lastBlock[[i]] <- c(strSizes[i, -N], strRows)
  maxStr <- do.call("pmax", lastBlock)
  for(i in 1:Q) lastBlock[[i]] <- maxStr[as.logical(lastBlock[[i]])]
  lastBlock <- c(lastBlock, list(X = strRows, y = strRows))
  list(N = N,                   # total number of rows in data
       ZXrows = N,              # no. of rows in array
       ZXcols = sum(ncols),     # no. of columns in array
       Q = Q,                   # no. of levels of random effects
       StrRows = strRows,       # no. of rows required for storage
       qvec = ncols * c(rep(1, Q), 0, 0), # lengths of random effects
                                        # no. of groups at each level
       ngrps = c(unlist(lapply(lastRow, length), N, N)),
                                        # offsets into DmHalf array by level
       DmOff = (c(0, cumsum(ncols^2)))[1:(Q+2)],
       ncol = ncols,            # no. of columns decomposed per level
                                        # no. of columns rotated per level
       nrot = (rev(c(0, cumsum(rev(ncols)))))[-1],
       ZXoff = offsets(N, ncols, lastRow), # offsets into ZXy
       ZXlen = lapply(lastRow, lengths), # lengths of ZXy groups
                                        # storage array offsets
       SToff = offsets(strRows, ncols, lastStr, triangle = TRUE),
                                        # decomposition offsets
       DecOff = offsets(strRows, ncols, lastBlock),
                                        # decomposition lengths
       DecLen = lapply(lastBlock, lengths)
       )
}

### Methods for standard generics

fitted.lme <-
  function(object, level = Q, asList = FALSE, ...)
{
#   Q <- object$dims$Q
#   val <- object[["fitted"]]
#   if (is.character(level)) {		# levels must be given consistently
#     nlevel <- match(level, names(val))
#     if (any(aux <- is.na(nlevel))) {
#       stop(paste("Nonexistent level(s)", level[aux]))
#     }
#     level <- nlevel
#   } else {				# assuming integers
#     level <- 1 + level
#   }
#   val2 <- napredict(object$na.action, val[, level])
#   if (length(level) == 1) {
#     grp.nm <- row.names(object[["groups"]])
#     grps <- as.character(object[["groups"]][, max(c(1, level - 1))])
#     if (asList) {
#       val <- as.list(split(val, ordered(grps, levels = unique(grps))))
#     } else {
#       grp.nm <- row.names(object[["groups"]])
#       val <- val2
#       names(val) <- grps[match(names(val), grp.nm)]
#     }
#     lab <- "Fitted values"
#     if (!is.null(aux <- attr(object, "units")$y)) {
#       lab <- paste(lab, aux)
#     }
#     attr(val, "label") <- lab
#     val
#   } else val2
}

formula.lme <- function(x, ...) eval(x$call$fixed)

fixef.lme <-
  function(object, ...) object$coefficients$fixed

newgetGroups.lme <-
  function(object, form, level = Q, data, sep)
{
#   Q <- object$dims$Q
#   val <- object[["groups"]][, level]
#   if (length(level) == 1) {		# single group
#     attr(val, "label") <- names(object[["groups"]])[level]
#   }
#   val
}

newgetGroupsFormula.lme <-
  function(object, asList = FALSE, sep)
{
#   getGroupsFormula(object$modelStruct$reStruct, asList)
}

getResponse.lme <-
  function(object, form)
{
#   val <- resid(object) + fitted(object)
#   if (is.null(lab <- attr(object, "labels")$y)) {
#     lab <- deparse(newgetResponseFormula(object)[[2]])
#   }
#   if (!is.null(aux <- attr(object, "units")$y)) {
#     lab <- paste(lab, aux)
#   }
#   attr(val, "label") <- lab
#   val
}

logLik.newlme <-
  function(object, REML, ...)
{
  p <- object$dims$ncol[object$dims$Q + 1]
  N <- object$dims$N
##  Np <- N - p
  estM <- object$method
  if (missing(REML)) REML <- estM == "REML"
  val <- object[["logLik"]]
  if (REML && (estM == "ML")) {			# have to correct logLik
    val <- val + (p * (log(2 * pi) + 1) + (N - p) * log(1 - p/N) +
		  sum(log(abs(svd(object$varFix)$d)))) / 2
  }
  if (!REML && (estM == "REML")) {	# have to correct logLik
    val <- val - (p * (log(2*pi) + 1) + N * log(1 - p/N) +
		  sum(log(abs(svd(object$varFix)$d)))) / 2
  }
  attr(val, "nall") <- N
  attr(val, "nobs") <- N - REML * p
  attr(val, "df") <- p + length(coef(object[["modelStruct"]])) + 1
  class(val) <- "logLik"
  val
}

print.newlme <-
  function(x, ...)
{
  dd <- x$dims
  if (inherits(x, "nlme")) {	# nlme object
    cat( "Nonlinear mixed-effects model fit by " )
    cat( ifelse( x$method == "REML", "REML\n", "maximum likelihood\n") )
    cat("  Model:", deparse(x$call$model),"\n")
  } else {				# lme objects
    cat( "Linear mixed-effects model fit by " )
    cat( ifelse( x$method == "REML", "REML\n", "maximum likelihood\n") )
  }
  cat("  Data:", deparse( x$call$data ), "\n")
  if (!is.null(x$call$subset)) {
    cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2]]),"\n")
  }
  cat("  Log-", ifelse(x$method == "REML", "restricted-", ""),
             "likelihood: ", format(x$logLik), "\n", sep = "")
  fixF <- x$call$fixed
  if (inherits(fixF, "formula") || is.call(fixF) || is.name(fixF)) {
    cat("  Fixed:", deparse(x$call$fixed), "\n")
  } else {
    cat("  Fixed:", deparse(lapply(fixF, function(el)
                                   as.name(deparse(el)))), "\n")
  }
  print(fixef(x))
  cat("\n")
  print(summary(x$modelStruct), sigma = x$sigma)
  cat("Number of Observations:", dd[["N"]])
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if ((lNgrps <- length(Ngrps)) == 1) {	# single nesting
    cat(Ngrps,"\n")
  } else {				# multiple nesting
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps),
			    c(lNgrps, lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps))
  }
  invisible(x)
}

print.summary.lme <-
  function(x, verbose = FALSE, ...)
{
  dd <- x$dims
  verbose <- verbose || attr(x, "verbose")
  if (inherits(x, "nlme")) {	# nlme object
    cat( "Nonlinear mixed-effects model fit by " )
    cat( ifelse( x$method == "REML", "REML\n", "maximum likelihood\n") )
    cat("  Model:", deparse(x$call$model),"\n")
  } else {				# lme objects
    cat( "Linear mixed-effects model fit by " )
    cat( ifelse( x$method == "REML", "REML\n", "maximum likelihood\n") )
  }
  method <- x$method
  cat(" Data:", deparse( x$call$data ), "\n")
  if (!is.null(x$call$subset)) {
    cat("  Subset:", deparse(asOneSidedFormula(x$call$subset)[[2]]),"\n")
  }
  print( data.frame( AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik),
                    row.names = " ") )
  if (verbose) { cat("Convergence at iteration:",x$numIter,"\n") }
  cat("\n")
  print(summary(x$modelStruct), sigma = x$sigma,
	reEstimates = x$coef$random, verbose = verbose)
  cat("Fixed effects: ")
  fixF <- x$call$fixed
  if (inherits(fixF, "formula") || is.call(fixF)) {
    cat(deparse(x$call$fixed), "\n")
  } else {
    cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))),
        "\n")
  }
  ## fixed effects t-table and correlations
  xtTab <- as.data.frame(x$tTable)
  wchPval <- match("p-value", names(xtTab))
  for(i in names(xtTab)[-wchPval]) {
    xtTab[, i] <- format(zapsmall(xtTab[, i]))
  }
  xtTab[,wchPval] <- format(round(xtTab[,wchPval], 4))
  if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0))) {
    levels(xtTab[, wchPval])[wchLv] <- "<.0001"
  }
  row.names(xtTab) <- dimnames(x$tTable)[[1]]
  print(xtTab)
  if (nrow(x$tTable) > 1) {
    corr <- x$corFixed
    class(corr) <- "correlation"
    print(corr,
	  title = " Correlation:",
	  ...)
  }
  cat("\nStandardized Within-Group Residuals:\n")
  print(x$residuals)
  cat("\nNumber of Observations:",x$dims[["N"]])
  cat("\nNumber of Groups: ")
  Ngrps <- dd$ngrps[1:dd$Q]
  if ((lNgrps <- length(Ngrps)) == 1) {	# single nesting
    cat(Ngrps,"\n")
  } else {				# multiple nesting
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps),
			    c(lNgrps, lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps))
  }
  invisible(x)
}

# qqnorm.lme <-
#   function(y, form = ~ resid(., type = "p"), abline = NULL,
#            id = NULL, idLabels = NULL, grid = FALSE, ...)
#   ## normal probability plots for residuals and random effects
# {
#   object <- y
#   if (!inherits(form, "formula")) {
#     stop("\"Form\" must be a formula")
#   }
#   ## constructing data
#   allV <- all.vars(newasOneFormula(form, id, idLabels))
#   allV <- allV[is.na(match(allV,c("T","F","TRUE","FALSE")))]
#   if (length(allV) > 0) {
#     data <- getData(object)
#     if (is.null(data)) {		# try to construct data
#       alist <- lapply(as.list(allV), as.name)
#       names(alist) <- allV
#       alist <- c(as.list(as.name("data.frame")), alist)
#       mode(alist) <- "call"
#       data <- eval(alist, sys.parent(1))
#     } else {
#       if (any(naV <- is.na(match(allV, names(data))))) {
# 	stop(paste(allV[naV], "not found in data"))
#       }
#     }
#   } else data <- NULL
#   ## argument list
#   dots <- list(...)
#   if (length(dots) > 0) args <- dots
#   else args <- list()
#   ## appending object to data
#   data <- as.list(c(as.list(data), . = list(object)))
# 
#   ## covariate - must always be present
#   covF <- newgetCovariateFormula(form)
#   .x <- eval(covF[[2]], data)
#   labs <- attr(.x, "label")
#   if (inherits(.x, "ranef.lme")) {      # random effects
#     type <- "reff"
#   } else {
#     if (!is.null(labs) && ((labs == "Standardized residuals") ||
#                            (labs == "Normalized residuals") ||
#                            (substring(labs, 1, 9) == "Residuals"))) {
#       type <- "res"                     # residuals
#     } else {
#       stop("Only residuals and random effects allowed")
#     }
#   }
#   if (is.null(args$xlab)) args$xlab <- labs
#   if (is.null(args$ylab)) args$ylab <- "Quantiles of standard normal"
#   if(type == "res") {			# residuals
#     fData <- qqnorm(.x, plot.it = FALSE)
#     data[[".y"]] <- fData$x
#     data[[".x"]] <- fData$y
#     dform <- ".y ~ .x"
#     if (!is.null(grp <- getGroupsFormula(form))) {
#       dform <- paste(dform, deparse(grp[[2]]), sep = "|")
#     }
#     if (!is.null(id)) {			# identify points in plot
#       id <-
#         switch(mode(id),
#                numeric = {
#                  if ((id <= 0) || (id >= 1)) {
#                    stop("Id must be between 0 and 1")
#                  }
#                  if (labs == "Normalized residuals") {
#                    as.logical(abs(resid(object, type="normalized"))
#                               > -qnorm(id / 2))
#                  } else {
#                    as.logical(abs(resid(object, type="pearson"))
#                               > -qnorm(id / 2))
#                  }
#                },
#                call = eval(asOneSidedFormula(id)[[2]], data),
#                stop("\"Id\" can only be a formula or numeric.")
#                )
#       if (is.null(idLabels)) {
#         idLabels <- getGroups(object)
#         if (length(idLabels) == 0) idLabels <- 1:object$dims$N
#         idLabels <- as.character(idLabels)
#       } else {
#         if (mode(idLabels) == "call") {
#           idLabels <-
#             as.character(eval(asOneSidedFormula(idLabels)[[2]], data))
#         } else if (is.vector(idLabels)) {
#           if (length(idLabels <- unlist(idLabels)) != length(id)) {
#             stop("\"IdLabels\" of incorrect length")
#           }
#           idLabels <- as.character(idLabels)
#         } else {
#           stop("\"IdLabels\" can only be a formula or a vector")
#         }
#       }
#     }
#   } else {				# random.effects
#     level <- attr(.x, "level")
#     std <- attr(.x, "standardized")
#     if (!is.null(effNams <- attr(.x, "effectNames"))) {
#       .x <- .x[, effNams, drop = FALSE]
#     }
#     nc <- ncol(.x)
#     nr <- nrow(.x)
#     fData <- lapply(as.data.frame(.x), qqnorm, plot.it = FALSE)
#     fData <- data.frame(.x = unlist(lapply(fData, function(x) x[["y"]])),
# 			.y = unlist(lapply(fData, function(x) x[["x"]])),
# 			.g = ordered(rep(names(fData),rep(nr, nc)),
#                         levels = names(fData)), check.names = FALSE)
#     dform <- ".y ~ .x | .g"
#     if (!is.null(grp <- getGroupsFormula(form))) {
#       dform <- paste(dform, deparse(grp[[2]]), sep = "*")
#       auxData <- data[is.na(match(names(data), "."))]
#     } else {
#       auxData <- list()
#     }
#     ## id and idLabels - need not be present
#     if (!is.null(id)) {			# identify points in plot
#       N <- object$dims$N
#       id <-
#         switch(mode(id),
#                numeric = {
#                  if ((id <= 0) || (id >= 1)) {
#                    stop("Id must be between 0 and 1")
#                  }
#                  aux <- ranef(object, level = level, standard = TRUE)
#                  as.logical(abs(c(unlist(aux))) > -qnorm(id / 2))
#                },
#                call = eval(asOneSidedFormula(id)[[2]], data),
#                stop("\"Id\" can only be a formula or numeric.")
#                )
#       if (length(id) == N) {
#         ## id as a formula evaluated in data
#         auxData[[".id"]] <- id
#       }
# 
#       if (is.null(idLabels)) {
#         idLabels <- rep(row.names(.x), nc)
#       } else {
#         if (mode(idLabels) == "call") {
#           idLabels <-
#             as.character(eval(asOneSidedFormula(idLabels)[[2]], data))
#         } else if (is.vector(idLabels)) {
#           if (length(idLabels <- unlist(idLabels)) != N) {
#             stop("\"IdLabels\" of incorrect length")
#           }
#           idLabels <- as.character(idLabels)
#         } else {
#           stop("\"IdLabels\" can only be a formula or a vector")
#         }
#       }
#       if (length(idLabels) == N) {
#         ## idLabels as a formula evaluated in data
#         auxData[[".Lid"]] <- idLabels
#       }
#     }
# 
#     if (length(auxData)) {		# need collapsing
#       auxData <- gsummary(as.data.frame(auxData),
#                           groups = getGroups(object, level = level))
#       auxData <- auxData[row.names(.x), , drop = FALSE]
# 
#       if (!is.null(auxData[[".id"]])) {
#         id <- rep(auxData[[".id"]], nc)
#       }
# 
#       if (!is.null(auxData[[".Lid"]])) {
#         idLabels <- rep(auxData[[".Lid"]], nc)
#       }
#       data <- cbind(fData, do.call("rbind", rep(list(auxData), nc)))
#     } else {
#       data <- fData
#     }
#   }
#   assign("id", if (is.null(id)) NULL else as.logical(as.character(id)))#,
#   #   where = 1)
#   assign("idLabels", as.character(idLabels))#, where = 1)
#   #assign("grid", grid, where = 1)
#   assign("abl", abline)#, where = 1)
#   if (is.null(args$strip)) {
#     args$strip <- function(...) strip.default(..., style = 1)
#   }
#   if (is.null(args$cex)) args$cex <- par("cex")
#   if (is.null(args$adj)) args$adj <- par("adj")
# 
#   args <- c(list(eval(parse(text = dform)),
#                  data = substitute(data)), args)
#   if (is.null(args$panel)) {
#     args <- c(list(panel = function(x, y, subscripts, ...){
#       x <- as.numeric(x)
#       y <- as.numeric(y)
#       dots <- list(...)
#       if (grid) panel.grid()
#       panel.xyplot(x, y, ...)
#       if (any(ids <- id[subscripts])){
#           ltext(x[ids], y[ids], idLabels[subscripts][ids],
#                 cex = dots$cex, adj = dots$adj)
#       }
#       if (!is.null(abl)) { if (length(abl) == 2) panel.abline(a = abl, ...) else panel.abline(h = abl, ...) }
#     }), args)
#   }
#   if(type == "reff" && !std) {
#     args[["scales"]] <- list(x = list(relation = "free"))
#   }
#   do.call("xyplot", as.list(args))
# }

# ranef.lme <-
#   ##  Extracts the random effects from an lme object.
#   ##  If aug.frame is true, the returned data frame is augmented with a
#   ##  values from the original data object, if available.  The variables
#   ##  in the original data are collapsed over the cluster variable by the
#   ##  function fun.
# function(object, augFrame = FALSE, level = 1:Q, data, which = 1:ncol(data),
# 	 FUN = mean, standard = FALSE , omitGroupingFactor = TRUE,
#          subset = NULL, ...)
# {
#   Q <- object$dims$Q
#   effects <- object$coefficients$random
#   if (Q > 1) {
#     grpNames <- t(array(rep(rev(names(effects)), Q), c(Q, Q)))
#     grpNames[lower.tri(grpNames)] <- ""
#     grpNames <-
#       rev(apply(grpNames, 1, function(x) paste(x[x != ""], collapse = " %in% ")))
#   } else {
#     grpNames <- names(effects)
#   }
#   effects <- effects[level]
#   grpNames <- grpNames[level]
#   if (standard) {
#     for (i in names(effects)) {
#       effects[[i]] <-
# 	t(t(effects[[i]]) / (object$sigma *
# 		     sqrt(diag(as.matrix(object$modelStruct$reStruct[[i]])))))
#     }
#   }
#   effects <- lapply(effects, as.data.frame)
#   if (augFrame) {
#     if (length(level) > 1) {
#       stop("Augmentation of random effects only available for single level")
#     }
#     effects <- effects[[1]]
#     effectNames <- names(effects)
#     if (missing(data)) {
#       data <- getData(object)
#     }
#     data <- as.data.frame(data)
#     if (is.null(subset)) {              # nlme case
#       subset <- eval(object$call[["naPattern"]])
#     } else {
#       subset <- asOneSidedFormula(as.list(match.call())[["subset"]])
#     }
#     if (!is.null(subset)) {
#       subset <- eval(subset[[2]], data)
#       data <- data[subset,  ,drop=FALSE]
#     }
#     data <- data[, which, drop = FALSE]
#     ## eliminating columns with same names as effects
#     data <- data[, is.na(match(names(data), effectNames)), drop = FALSE]
#     grps <- as.character(object[["groups"]][, level])
#     data <- gsummary(data, FUN = FUN, groups = grps)
#     if (omitGroupingFactor) {
#       data <-
# 	data[, is.na(match(names(data), names(object$modelStruct$reStruct))),
# 	      drop = FALSE]
#     }
#     if (length(data) > 0) {
#       effects <- cbind(effects, data[row.names(effects),, drop = FALSE])
#     }
#     attr(effects, "effectNames") <- effectNames
#   } else {
#     effects <- lapply(effects,
#                       function(el) {
#                         attr(el, "effectNames") <- names(el)
#                         el
#                       })
#     if (length(level) == 1) effects <- effects[[1]]
#   }
#   attr(effects, "label") <-
#     if (standard) {
#       "Standardized random effects"
#     } else {
#       "Random effects"
#     }
#   attr(effects, "level") <- max(level)
#   attr(effects, "standardized") <- standard
#   attr(effects, "grpNames") <- grpNames
#   class(effects) <- c("ranef.lme", class(effects))
#   effects
# }

# residuals.lme <-
#   function(object, level = Q, type = c("response", "pearson", "normalized"),
#            asList = FALSE, ...)
# 
# {
#   type <- match.arg(type)
#   Q <- object$dims$Q
#   val <- object[["residuals"]]
#   if (is.character(level)) {		# levels must be given consistently
#     nlevel <- match(level, names(val))
#     if (any(aux <- is.na(nlevel))) {
#       stop(paste("Nonexistent level(s)", level[aux]))
#     }
#     level <- nlevel
#   } else {				# assuming integers
#     level <- 1 + level
#   }
#   if (type != "response") {		# standardize
#     ## have to standardize properly for when corStruct neq NULL
#     val <- val[, level]/attr(val, "std")
#   } else {
#     val <- val[, level]
#   }
#   if (type == "normalized") {
#     if (!is.null(cSt <- object$modelStruct$corStruct)) {
#       ## normalize according to inv-trans factor
#       val <- recalc(cSt, list(Xy = as.matrix(val)))$Xy[, 1:length(level)]
#     } else {                            # will just standardized
#       type <- "pearson"
#     }
#   }
#   if (length(level) == 1) {
#     grps <- as.character(object[["groups"]][, max(c(1, level - 1))])
#     if (asList) {
#       val <- as.list(split(val, ordered(grps, levels = unique(grps))))
#     } else {
#       grp.nm <- row.names(object[["groups"]])
#       val <-naresid(object$na.action, val)
#       names(val) <- grps[match(names(val), grp.nm)]
#     }
#     attr(val, "label") <-
#       switch(type,
#              response = {
#                lab <- "Residuals"
#                if (!is.null(aux <- attr(object, "units")$y)) {
#                  lab <- paste(lab, aux)
#                }
#                lab
#              },
#              pearson = "Standardized residuals",
#              normalized = "Normalized residuals"
#              )
#     val
#   } else naresid(object$na.action, val)
# }

summary.newlme <- function(object, adjustSigma = TRUE, verbose = FALSE, ...)
{
  ##  variance-covariance estimates for fixed effects
  fixed <- fixef(object)
  stdFixed <- sqrt(diag(as.matrix(object$varFix)))
  object$corFixed <- array(t(object$varFix/stdFixed)/stdFixed,
                           dim(object$varFix), list(names(fixed),names(fixed)))
  if (object$method == "ML" && adjustSigma == TRUE) {
    stdFixed <-
      sqrt(object$dims$N/(object$dims$N - length(stdFixed))) * stdFixed
  }
  ## fixed effects coefficients, std. deviations and t-ratios
  ##
  tTable <- data.frame(fixed, stdFixed, object$fixDF[["X"]],
                       fixed/stdFixed, fixed)
  dimnames(tTable)<-
    list(names(fixed),c("Value", "Std.Error", "DF", "t-value", "p-value"))
  tTable[, "p-value"] <- 2 * pt(-abs(tTable[,"t-value"]), tTable[,"DF"])
  object$tTable <- as.matrix(tTable)
  ##
  ## residuals
  ##
  resd <- resid(object, type = "pearson")
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE) # might have NAs from na.exclude
    names(resd) <- c("Min","Q1","Med","Q3","Max")
  }
  object$residuals <- resd
  ##
  ## generating the final object
  ##
  aux <- logLik(object)
  object$BIC <- BIC(aux)
  object$AIC <- AIC(aux)
  attr(object, "oClass") <- class(object)
  attr(object, "verbose") <- verbose
  class(object) <- c("summary.lme", class(object))

  object
}

# based on R's update.default
update.newlme <-
    function (object, fixed., ..., evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
	stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(fixed.))
	call$fixed <- update.formula(formula(object), fixed.)
    if(length(extras) > 0) {
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

#update.lme <-
#  function(object, fixed, data, random, correlation, weights, subset,
#           method, na.action, control, contrasts, ...)
#{
#  thisCall <- as.list(match.call())[-(1:2)]
#  if (is.null(nextCall <- object$origCall) ||
#      !is.null(thisCall$fixed) ||
#      is.null(thisCall$random)) {
#    nextCall <- object$call
#  }
#  nextCall <- as.list(nextCall)[-1]
#  if (is.null(thisCall$random)  && is.null(thisCall$subset)) {
#    ## no changes in ranef model and no subsetting
#    thisCall$random <- object$modelStruct$reStruct
#  }
#  if (is.na(match("correlation", names(thisCall))) &&
#      !is.null(thCor <- object$modelStruct$corStruct)) {
#    thisCall$correlation <- thCor
#  }
#  if (is.na(match("weights", names(thisCall))) &&
#      !is.null(thWgt <- object$modelStruct$varStruct)) {
#    thisCall$weights <- thWgt
#  }
#    argNams <- unique( c(names(nextCall), names(thisCall)) )
#    args <- vector("list", length(argNams))
#    names(args) <- argNams
#    args[ names(nextCall) ] <- nextCall
#    nextCall <- args
#  if (!is.null(thisCall$fixed)) {
#    thisCall$fixed <- update(as.formula(nextCall$fixed), fixed)
#  }
#  nextCall[names(thisCall)] <- thisCall
#  do.call("lme", nextCall)
#}

# Variogram.lme <-
#   function(object, distance, form = ~1,
#            resType = c("pearson", "response", "normalized"),
#            data, na.action = na.fail, maxDist, length.out = 50,
#            collapse = c("quantiles", "fixed", "none"), nint = 20, breaks,
#            robust = FALSE, metric = c("euclidean", "maximum", "manhattan"),
#            ...)
# {
#   resType <- match.arg(resType)
#   Q <- object$dims$Q
#   grps <- getGroups(object, level = Q)
#   ## checking if object has a corSpatial element
#   csT <- object$modelStruct$corStruct
#   wchRows <- NULL
#   if (missing(distance)) {
#     if (missing(form) && inherits(csT, "corSpatial")) {
#       distance <- newgetCovariate(csT)
#     } else {
#       metric <- match.arg(metric)
#       if (missing(data)) {
#         data <- getData(object)
#       }
#       if (is.null(data)) {			# will try to construct
#         allV <- all.vars(form)
#         if (length(allV) > 0) {
#           alist <- lapply(as.list(allV), as.name)
#           names(alist) <- allV
#           alist <- c(as.list(as.name("data.frame")), alist)
#           mode(alist) <- "call"
#           data <- eval(alist, sys.parent(1))
#         }
#       }
#       covForm <- newgetCovariateFormula(form)
#       if (length(all.vars(covForm)) > 0) {
#         if (attr(terms(covForm), "intercept") == 1) {
#           covForm <-
#             eval(parse(text = paste("~", c_deparse(covForm[[2]]),"-1",sep="")))
#         }
#         covar <- model.frame(covForm, data, na.action = na.action)
#         ## making sure grps is consistent
#         wchRows <- !is.na(match(row.names(data), row.names(covar)))
#         grps <- grps[wchRows, drop = TRUE]
#         covar <- as.data.frame(unclass(model.matrix(covForm, covar)))
#       } else {
#         covar <-
#           data.frame(dist = unlist(tapply(rep(1, nrow(data)), grps, cumsum)))
#       }
#       covar <- split(covar, grps)
#       ## getting rid of 1-observation groups
#       covar <- covar[sapply(covar, function(el) nrow(as.matrix(el))) > 1]
#       distance <- lapply(covar,
#                          function(el, metric) dist(as.matrix(el), metric),
#                          metric = metric)
#     }
#   }
#   res <- resid(object, type = resType)
#   if (!is.null(wchRows)) {
#     res <- res[wchRows]
#   }
#   res <- split(res, grps)
#   res <- res[sapply(res, length) > 1] # no 1-observation groups
#   levGrps <- levels(grps)
#   val <- vector("list", length(levGrps))
#   names(val) <- levGrps
#   for(i in levGrps) {
#     val[[i]] <- Variogram(res[[i]], distance[[i]])
#   }
#   val <- do.call("rbind", val)
#   if (!missing(maxDist)) {
#     val <- val[val$dist <= maxDist, ]
#   }
#   collapse <- match.arg(collapse)
#   if (collapse != "none") {             # will collapse values
#     dst <- val$dist
#     udist <- sort(unique(dst))
#     ludist <- length(udist)
#     if (!missing(breaks)) {
#       if (min(breaks) > udist[1]) {
#         breaks <- c(udist[1], breaks)
#       }
#       if (max(breaks) < udist[2]) {
#         breaks <- c(breaks, udist[2])
#       }
#       if (!missing(nint) && nint != (length(breaks) - 1)) {
#         stop("'nint' is not consistent with 'breaks'")
#       }
#       nint <- length(breaks) - 1
#     }
#     if (nint < ludist) {
#       if (missing(breaks)) {
#         if (collapse == "quantiles") {    # break into equal groups
#           breaks <- unique(quantile(dst, seq(0, 1, 1/nint)))
#         } else {                          # fixed length intervals
#           breaks <- seq(udist[1], udist[length(udist)], length = nint + 1)
#         }
#       }
#       cutDist <- cut(dst, breaks)
#     } else {
#       cutDist <- dst
#     }
#     val <- lapply(split(val, cutDist),
#                   function(el, robust) {
#                     nh <- nrow(el)
#                     vrg <- el$variog
#                     if (robust) {
#                       vrg <- ((mean(vrg^0.25))^4)/(0.457+0.494/nh)
#                     } else {
#                       vrg <- mean(vrg)
#                     }
#                     dst <- median(el$dist)
#                     data.frame(variog = vrg, dist = dst)
#                   }, robust = robust)
#     val <- do.call("rbind", val)
#     val$n.pairs <- as.vector(table(na.omit(cutDist)))
#     val <- na.omit(val)                 # getting rid of NAs
#   }
#   row.names(val) <- 1:nrow(val)
#   if (inherits(csT, "corSpatial") && resType != "normalized") {
#     ## will keep model variogram
#     if (resType == "pearson") {
#       sig2 <- 1
#     } else {
#       sig2 <- object$sigma^2
#     }
#     attr(val, "modelVariog") <-
#       Variogram(csT, sig2 = sig2, length.out = length.out)
#   }
#   attr(val, "collapse") <- collapse != "none"
#   class(val) <- c("Variogram", "data.frame")
#   val
# }

###*### newlmeStruct - a model structure for lme fits

newlmeStruct <-
  ## constructor for newlmeStruct objects
  function(reStruct, corStruct = NULL, varStruct = NULL)
{
  val <- list(reStruct = reStruct, corStruct = corStruct,
              varStruct = varStruct)
  val <- val[!sapply(val, is.null)]	# removing NULL components
  attr(val, "settings") <- attr(val$reStruct, "settings")
  class(val) <- c("newlmeStruct", "modelStruct")
  
  val
}

##*## newlmeStruct methods for standard generics

fitted.newlmeStruct <-
  function(object, level = Q, conLin = attr(object, "conLin"),
           lmeFit = attr(object, "lmeFit"), ...)
{
  if (is.null(conLin)) {
    stop("No condensed linear model")
  }
  if (is.null(lmeFit)) {
    stop("No fitted lme object")
  }
  dd <- conLin$dims
  Q <- dd$Q
  Qp1 <- Q + 1
  nc <- dd$ncol
  fit <- array(0, c(dd$N, Qp1),
       list(dimnames(conLin$Xy)[[1]], c("fixed", rev(names(object$reStruct)))))
  ZXstart <- rev(cumsum(c(1, nc[1:Q])))
  ZXend <- rev(cumsum(nc[1:Qp1]))
  ZXlen <- dd$ZXlen[Q:1]
  ZXngrps <- dd$ngrps[Q:1]
  ZXb <- lmeFit$b
  nc <- nc[Q:1]

  fit[, "fixed"] <-			# population fitted values
    conLin$Xy[, ZXstart[1]: ZXend[1], drop = FALSE] %*% lmeFit$beta

  for(i in 1:Q) {
    j <- i + 1
    fit[, j] <- fit[, i] +
      (conLin$Xy[, ZXstart[j]:ZXend[j], drop = FALSE] *
       ZXb[[i]][rep(1:ZXngrps[i], ZXlen[[i]]),,drop = FALSE]) %*% rep(1, nc[i])
  }
  ## this is documented to return a vector for one level, matrix for more.
  ## So it should be a matrix if there is only one row, but not if
  ## there is only one columns.
  if(length(level) > 1) fit[, level + 1, drop = FALSE] else fit[, level+1]
}

newInitialize.newlmeStruct <-
  function(object, data, groups, conLin = attr(object, "conLin"),
	   control= list(niterEM = 20, gradHess = TRUE), ...)
{
  object[] <- lapply(object, newInitialize, data, conLin, control)
  theta <- lapply(object, coef)
  len <- unlist(lapply(theta, length))
  num <- seq_along(len)
  if (sum(len) > 0) {
    pmap <- outer(rep(num, len), num, "==")
  } else {
    pmap <- array(FALSE, c(1, length(len)))
  }
  dimnames(pmap) <- list(NULL, names(object))
  attr(object, "pmap") <- pmap
  if (length(object) == 1  &&           # only reStruct
      all(attr(object, "settings")[-(1:3)] >= 0) && # known pdMat class
      control[["gradHess"]]) {
    ## can use numerical derivatives
    attr(object, "settings")[2:3] <- c(0, 1)
    class(object) <- c("newlmeStructInt", class(object))
  }
  if (needUpdate(object)) {
    attr(object, "lmeFit") <- newMEestimate(object, groups)
    update(object, data)
  } else {
    object
  }
}

logLik.newlmeStruct <-
  function(object, Pars, conLin = attr(object, "conLin"), ...)
{
  coef(object) <- Pars			# updating parameter values
  newrecalc(object, conLin)[["logLik"]]	# updating conLin
}

logLik.newlmeStructInt <-
  function(object, Pars, conLin = attr(object, "conLin"), ...)
{
  ## logLik for objects with reStruct parameters only, with
  ## internally defined class
  q <- length(Pars)
  aux <- .C("mixed_loglik",
	    as.double(conLin[["Xy"]]),
	    as.integer(unlist(conLin$dims)),
	    as.double(Pars),
	    as.integer(attr(object, "settings")),
	    val = double(1 + q * (q + 1)),
	    double(1))[["val"]]
  val <- aux[1]
  attr(val, "gradient") <- -aux[1 + (1:q)]
  attr(val, "hessian") <- -array(aux[-(1:(q+1))], c(q, q))
  val
}

residuals.newlmeStruct <-
  function(object, level = Q, conLin = attr(object, "conLin"),
           lmeFit = attr(object, "lmeFit"), ...)
{
  Q <- conLin$dims$Q
  conLin$Xy[, conLin$dims$ZXcols] - fitted(object, level, conLin, lmeFit)
}

varWeights.newlmeStruct <-
  function(object)
{
  if (is.null(object$varStruct)) rep(1, attr(object, "conLin")$dims$N)
  else varWeights(object$varStruct)
}

## Auxiliary control functions

newlmeScale <- function(start)
#
# function used to set the scale inside ms(), for lme() and nlme()
# calls
#
{
  scale <- abs(start)
  nonzero <- scale > 0
  if (any(nonzero)) {
    scale[nonzero] <- 1/scale[nonzero]
    scale[!nonzero] <- median(scale[nonzero])
  }
  else {
    scale <- rep(1, length(scale))
  }
  scale
}

newlmeControl <-
  ## Control parameters for lme
  function(maxIter = 50, msMaxIter = 50, tolerance = 1e-6, niterEM = 25, msMaxEval = 200,
	   msTol = 1e-7, msScale = lmeScale, msVerbose = FALSE,
           returnObject = FALSE, gradHess = TRUE, apVar = TRUE,
	   .relStep = (.Machine$double.eps)^(1/3), minAbsParApVar = 0.05,
           nlmStepMax = 100.0, opt = c("nlminb", "optim"),
	   optimMethod = "BFGS", natural = TRUE)
{
  list(maxIter = maxIter, msMaxIter = msMaxIter, tolerance = tolerance,
       niterEM = niterEM, msMaxEval = msMaxEval, msTol = msTol, msScale = msScale,
       msVerbose = msVerbose, returnObject = returnObject,
       gradHess = gradHess , apVar = apVar, .relStep = .relStep,
       nlmStepMax = nlmStepMax, opt = match.arg(opt),
       optimMethod = optimMethod,
       minAbsParApVar = minAbsParApVar, natural = natural)
}


### Constructor
### There is no constructor function for this class (i.e. no function
### called modelStruct) because the class is virtual.
### Objects inheriting from this class are required to have a "conLin"
### (condensed linear model) attribute and a "pmap" (parameter map)
### attribute

###*# Methods for standard generics

coef.modelStruct <-
  function(object, unconstrained = TRUE, ...)
{
  unlist(lapply(object, coef, unconstrained))
}

"coef<-.modelStruct" <-
  function(object, ..., value)
{
  value <- as.numeric(value) 
  parMap <- attr(object, "pmap")
  for(i in names(object)) {
    if (any(parMap[,i])) {
      coef(object[[i]]) <- value[parMap[,i]]
    }
  }
  object
}

formula.modelStruct <-
  function(x, ...)
{
  lapply(x, formula)
}

needUpdate.modelStruct <-
  function(object) any(unlist(lapply(object, needUpdate)))

print.modelStruct <-
  function(x, ...)
{
  for(i in names(x)) {
    if ((length(aux <- coef(x[[i]]))) > 0) {
      cat(paste(i, " parameters:\n"))
      print(aux)
    }
  }
  invisible(x)
}

print.summary.modelStruct <-
  function(x, ...)
{
  lapply(x, print, ...)
  invisible(x)
}

newrecalc.modelStruct <-
  function(object, conLin = attr(object, "conLin"), ...)
{
  for(i in rev(seq_along(object))) {
    conLin <- newrecalc(object[[i]], conLin)
    NULL
  }
  conLin
}

summary.modelStruct <-
  function(object, ...)
{
  value <- lapply(object, summary)
  class(value) <- "summary.modelStruct"
  value
}

## will not work as it is. fitted needs more than one argument!
update.modelStruct <-
  function(object, data, ...)
{
  if (needUpdate(object)) {
    object[] <- lapply(object, update, c(list("." = object), as.list(data)))
  }
  object
}

## newFunc.R

allCoef <-
  ## Combines different coefficient vectors into one vector, keeping track
  ## of which coefficients came from which object
  function(..., extract = coef)
{
  dots <- list(...)
  theta <- lapply(dots, extract)
  len <- unlist(lapply(theta, length))
  num <- seq_along(len)
  if (sum(len) > 0) {
    which <- outer(rep(num, len), num, "==")
  } else {
    which <- array(FALSE, c(1, length(len)))
  }
  cnames <- unlist(as.list(sys.call()[-1]))
  dimnames(which) <- list(NULL, cnames[cnames != substitute(extract)])
  theta <- unlist(theta)
  attr(theta, "map") <- which
  theta
}

newallVarsRec <-
  ## Recursive version of all.vars
  function(object)
{
  if (is.list(object)) {
    unlist(lapply(object, newallVarsRec))
  } else {
    all.vars(object)
  }
}

newasOneFormula <-
  ## Constructs a linear formula with all the variables used in a
  ## list of formulas, except for the names in omit
  function(..., omit = c(".", "pi"))
{
  names <- unique(newallVarsRec((list(...))))
  names <- names[is.na(match(names, omit))]
  if (length(names)) {
    eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  } else NULL
}

 compareFits <-
   ## compares coeffificients from different fitted objects
   function(object1, object2, which = 1:ncol(object1))
 {
#   dn1 <- dimnames(object1)
#   dn2 <- dimnames(object2)
#   aux <- rep(NA, length(dn1[[1]]))
#   if (any(aux1 <- is.na(match(dn2[[2]], dn1[[2]])))) {
#     object1[,dn2[[2]][aux1]] <- aux
#   }
#   if (any(aux1 <- is.na(match(dn1[[2]], dn2[[2]])))) {
#     object2[,dn1[[2]][aux1]] <- aux
#   }
#   dn1 <- dimnames(object1)
#   c1 <- deparse(substitute(object1))
#   c2 <- deparse(substitute(object2))
#   if (any(sort(dn1[[1]]) != sort(dn2[[1]]))) {
#     stop("Objects must have coefficients with same row names")
#   }
#   ## putting object2 in same order
#   object2 <- object2[dn1[[1]], dn1[[2]], drop = FALSE]
#   object1 <- object1[, which, drop = FALSE]
#   object2 <- object2[, which, drop = FALSE]
#   dn1 <- dimnames(object1)
#   dm1 <- dim(object1)
#   out <- array(0, c(dm1[1], 2, dm1[2]), list(dn1[[1]], c(c1,c2), dn1[[2]]))
#   for(i in dn1[[2]]) {
#     out[,,i] <- cbind(object1[[i]], object2[[i]])
#   }
#   class(out) <- "compareFits"
#   out
 }

newfdHess <- function(pars, fun, ..., .relStep = (.Machine$double.eps)^(1/3),
                    minAbsPar = 0)
   ## Use a Koschal design to establish a second order model for the response
 {
  pars <- as.numeric(pars)
  npar <- length(pars)
  incr <- ifelse( abs(pars) <= minAbsPar, minAbsPar * .relStep,
                 abs(pars) * .relStep )
  baseInd <- diag(npar)
  frac <- c(1, incr, incr^2)
  cols <- list(0, baseInd, -baseInd)
  for ( i in seq_along(pars)[ -npar ] ) {
    cols <- c( cols, list( baseInd[ , i ] + baseInd[ , -(1:i) ] ) )
    frac <- c( frac, incr[ i ] * incr[ -(1:i) ] )
  }
  indMat <- do.call( "cbind", cols)
  shifted <- pars + incr * indMat
  indMat <- t(indMat)
  Xcols <- list(1, indMat, indMat^2)
  for ( i in seq_along(pars)[ - npar ] ) {
    Xcols <- c( Xcols, list( indMat[ , i ] * indMat[ , -(1:i) ] ) )
  }
  coefs <- solve( do.call( "cbind", Xcols ) , apply(shifted, 2, fun, ...) )/frac
  Hess <- diag( coefs[ 1 + npar + seq_along(pars) ], ncol = npar )
  Hess[ row( Hess ) > col ( Hess ) ] <- coefs[ -(1:(1 + 2 * npar)) ]
  list( mean = coefs[ 1 ], gradient = coefs[ 1 + seq_along(pars) ],
       Hessian = ( Hess + t(Hess) ) )
}

gapply <-
  ## Apply a function to the subframes of a data.frame
  ## If "apply" were generic, this would be the method for groupedData
  function(object, which, FUN, form = formula(object), level,
           groups = newgetGroups(object, form, level), ...)
{
  if (!inherits(object, "data.frame")) {
    stop("Object must inherit from data.frame")
  }
  ## Apply a function to the subframes of a groupedData object
  if (missing(groups)) {                # formula and level are required
    if (!inherits(form, "formula")) {
      stop("\"Form\" must be a formula")
    }
    if (is.null(grpForm <- newgetGroupsFormula(form, asList = TRUE))) {
      ## will use right hand side of form as groups formula
      grpForm <- newsplitFormula(asOneSidedFormula(form[[length(form)]]))
    }
    if (missing(level)) level <- length(grpForm)
    else if (length(level) != 1) {
      stop("Only one level allowed in gapply")
    }
    groups <- groups                    # forcing evaluation
  }
  if (!missing(which)) {
    switch(mode(which),
           character = {
             wchNot <- is.na(match(which, names(object)))
             if (any(wchNot)) {
               stop(paste(paste(which[wchNot], collapse = ","),
                          "not matched"))
             }
           },
           numeric = {
             if (any(is.na(match(which, 1:ncol(object))))) {
               stop("Which must be between 1 and", ncol(object))
             }
           },
           stop("Which can only be character or integer.")
           )
    object <- object[, which, drop = FALSE]
  }
  val <- lapply(split(object, groups), FUN, ...)
  if (is.atomic(val[[1]]) && length(val[[1]]) == 1) {
    val <- unlist(val)
  }
  val
}

newgetCovariateFormula <-
  function(object)
{
  ## Return the primary covariate formula as a one sided formula
  form <- formula(object)
  if (!(inherits(form, "formula"))) {
    stop("formula(object) must return a formula")
  }
  form <- form[[length(form)]]
  if (length(form) == 3 && form[[1]] == as.name("|")){ # conditional expression
    form <- form[[2]]
  }
  eval(substitute(~form))
}

newgetResponseFormula <-
  function(object)
{
  ## Return the response formula as a one sided formula
  form <- formula(object)
  if (!(inherits(form, "formula") && (length(form) == 3))) {
    stop("\"Form\" must be a two sided formula")
  }
  eval(parse(text = paste("~", deparse(form[[2]]))))
}

gsummary <-
  ## Summarize an object according to the levels of a grouping factor
  ##
  function(object, FUN = function(x) mean(x, na.rm = TRUE),
           omitGroupingFactor = FALSE,
	   form = formula(object), level,
	   groups = newgetGroups(object, form , level),
	   invariantsOnly = FALSE, ...)
{
  if (!inherits(object, "data.frame")) {
    stop("Object must inherit from data.frame")
  }
  if (missing(groups)) {                # formula and level are required
    if (!inherits(form, "formula")) {
      stop("\"Form\" must be a formula")
    }
    if (is.null(grpForm <- newgetGroupsFormula(form, asList = TRUE))) {
      ## will use right hand side of form as groups formula
      grpForm <- newsplitFormula(asOneSidedFormula(form[[length(form)]]))
    }
    if (missing(level)) level <- length(grpForm)
    else if (length(level) != 1) {
      stop("Only one level allowed in gsummary")
    }
  }
  gunique <- unique(groups)
  firstInGroup <- match(gunique, groups)
  asFirst <- firstInGroup[match(groups, gunique)]
  value <- as.data.frame(object[firstInGroup, , drop = FALSE])
  row.names(value) <- as.character(gunique)
  value <- value[as.character(sort(gunique)), , drop = FALSE]
  varying <- unlist(lapply(object,
			   function(column, frst) {
			     aux <- as.character(column)
			     any(!identical(aux, aux[frst]))
			   },
			   frst = asFirst))
  if (any(varying) && (!invariantsOnly)) { # varying wanted
    Mode <- function(x) {
      aux <- table(x)
      names(aux)[match(max(aux), aux)]
    }
    if (data.class(FUN) == "function") {	# single function given
      FUN <- list(numeric = FUN, ordered = Mode, factor = Mode)
    } else {
      if (!(is.list(FUN) &&
	   all(sapply(FUN, data.class) == "function"))) {
	stop("FUN can only be a function or a list of functions")
      }
      auxFUN <- list(numeric = mean, ordered = Mode, factor = Mode)
      aux <- names(auxFUN)[is.na(match(names(auxFUN), names(FUN)))]
      if (length(aux) > 0) FUN[aux] <- auxFUN[aux]
    }
    for(nm in names(object)[varying]) {
      ## dClass <- data.class(object[[nm]])
      ## The problem here is that dclass may find an irrelevant class,
      ## e.g. Hmisc's "labelled"
      dClass <- if(is.ordered(object[[nm]])) "ordered" else
	        if(is.factor(object[[nm]])) "factor" else mode(object[[nm]])
      if (dClass == "numeric") {
	value[[nm]] <- as.vector(tapply(object[[nm]], groups, FUN[["numeric"]],...))
      } else {
	value[[nm]] <-
	  as.vector(tapply(as.character(object[[nm]]), groups, FUN[[dClass]]))
        if (inherits(object[,nm], "ordered")) {
          value[[nm]] <- ordered(value[,nm], levels = levels(object[,nm]))[drop = TRUE]
        } else {
          value[[nm]] <- factor(value[,nm], levels = levels(object[,nm]))[drop = TRUE]
        }
      }
    }
  } else {				# invariants only
    value <- value[, !varying, drop = FALSE]
  }
  if (omitGroupingFactor) {
    if (is.null(form)) {
      stop("Cannot omit grouping factor without \"form\"")
    }
    grpForm <- newgetGroupsFormula(form, asList = TRUE)
    if (missing(level)) level <- length(grpForm)
    grpNames <- names(grpForm)[level]
    whichKeep <- is.na(match(names(value), grpNames))
    if (any(whichKeep)) {
      value <- value[ , whichKeep, drop = FALSE]
    } else {
      return(NULL);
    }
  }
  value
}

newsplitFormula <-
  ## split, on the nm call, the rhs of a formula into a list of subformulas
  function(form, sep = "/")
{
  if (inherits(form, "formula") ||
      mode(form) == "call" && form[[1]] == as.name("~"))
    return(newsplitFormula(form[[length(form)]], sep = sep))
  if (mode(form) == "call" && form[[1]] == as.name(sep))
    return(do.call("c", lapply(as.list(form[-1]), newsplitFormula, sep = sep)))
  if (mode(form) == "(") return(newsplitFormula(form[[2]], sep = sep))
  if (length(form) < 1) return(NULL)
  list(asOneSidedFormula(form))
}

###    New generics used with corStruct, varFunc, groupedData, and reStruct


newACF <-
  ## autocorrelation function - needed not exist if acf were generic
  function(object, maxLag, ...) UseMethod("newACF")

newBIC <-
  ## Return the object's value of the Bayesian Information Criterion
  function(object, ...) UseMethod("newBIC")

newasTable <-
  ## Return the object in a tabular form
  function(object) UseMethod("newasTable")

newaugPred <-
  ## Return the data used to fit the model augmented with the predictions
  function(object, primary = NULL, minimum = min(primary),
           maximum = max(primary), length.out = 51, ...) UseMethod("newaugPred")

"coef<-" <- "coefficients<-" <-
  ## Assignment of the unconstrained parameter
  function(object, ..., value) UseMethod("coef<-")

newcollapse <-
  ## collapse a data frame according to a factor, or several nested factors
  function(object, ...) UseMethod("newcollapse")

newcomparePred <-
  ## compare predictions from different fitted objects
  function(object1, object2, primary = NULL,
	   minimum = min(primary), maximum = max(primary),
	   length.out = 51, level = NULL, ...) UseMethod("newcomparePred")

"covariate<-" <-
  ## Assignment of the primary covariate
  function(object, value) UseMethod("covariate<-")

newDim <-
  ## Extract dimensions of an object. Not needed if "dims" were generic
  function(object, ...) UseMethod("newDim")

newfixed.effects <-
  ## Generic extractor for estimates of fixed effects
  function(object, ...) UseMethod("newfixef")

newfixef <-
  ## Short form for generic extractor for estimates of fixed effects
  function(object, ...) UseMethod("newfixef")

newgetCovariate <-
  ## Return the primary covariate associated with object according to form
  function(object, form = formula(object), data)
    UseMethod("newgetCovariate")

newgetData <-
  ## Return the data.frame used to fit an object, if any was given in
  ## the call that produced it
  function(object) UseMethod("newgetData")

newgetGroups <-
  ## Return the groups associated with object according to form.
  function(object, form = formula(object), level, data, sep = "/")
    UseMethod("newgetGroups")

newgetGroupsFormula <-
  ## Return the formula(s) for the groups associated with object.
  ## The result is a one-sided formula unless asList is TRUE in which case
  ## it is a list of formulas, one for each level.
  function(object, asList = FALSE, sep = "/")
    UseMethod("newgetGroupsFormula")

getResponse <-
  ## Return the response associated with object according to form.
  function(object, form = formula(object))
    UseMethod("getResponse")

newisBalanced <-
  ## Check for balance, especially in a groupedData object
  function(object, countOnly = FALSE, level) UseMethod("newisBalanced")

newisInitialized <-
  ## Determine if the object has been assigned a value
  function(object) UseMethod("newisInitialized")

newInitialize <-
  ## newInitialize  objects
  function(object, data, ...) UseMethod("newInitialize")

newintervals <-
  ## generate confidence intervals for the parameters in object
  function(object, level = 0.95, ...) UseMethod("newintervals")

newlogDet <-
  ## Returns the negative of the sum of the logarithm of the determinant
  function(object, ...) UseMethod("newlogDet")

"matrix<-" <-
  ## Assignment of the matrix in an object representing special types of matrices
  function(object, value) UseMethod("matrix<-")

Names <-
  ## Extract names of an object. Not needed if "names" were generic
  function(object, ...) UseMethod("Names")

"Names<-" <-
  ## Assignment of names. Not needed if "names<-" were generic
  function(object, ..., value) UseMethod("Names<-")

newneedUpdate <-
  ## Checks if model plug-in needs to be updated after an estimation cycle
  function(object) UseMethod("newneedUpdate")

#pruneLevels <-
#  ## Returns the factor with the levels attribute truncated to only those
#  ## levels occuring in the factor
#  function(object) UseMethod("pruneLevels")

newrandom.effects <-
  ## Generic function for extracting the random effects
  ## If aug.frame is true, the returned data frame is augmented with
  ## values from the original data object, if available.  The variables
  ## in the original data are collapsed over the groups variable by the
  ## function fun.
  function(object, ...) UseMethod("newranef")

newranef <-
  ## Short form for generic function for extracting the random effects
  function(object, ...) UseMethod("newranef")

newrecalc <-
  ## Recalculate condensed linear object, according to model plug-in
  function(object, conLin, ...) UseMethod("newrecalc")

newVariogram <-
  ## calculates variogram of a vector according to a distance matrix
  function(object, distance, ...)
  UseMethod("newVariogram")


###      Methods for generics from newGenerics.q for some standard classes


##*## Methods for some of the generics in newGenerics.q for standard classes

newBIC.logLik <-
  ## BIC for logLik objects
  function(object, ...)
{
  -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
}

newBIC.lm <- 
 ## BIC for various fitted objects
  function(object, ...)
{
#   if(nargs() > 1) {
#     object <- list(object, ...)
#     val <- lapply(object, logLik)
#     val <-
#       as.data.frame(t(sapply(val, function(el) c(attr(el, "df"), BIC(el)))))
#     names(val) <- c("df", "BIC")
#     row.names(val) <- as.character(match.call()[-1])
#     return(val)
#   }
#   BIC(logLik(object))
}

newDim.default <- function(object, ...) dim(object)

newgetCovariate.data.frame <-
  function(object, form = formula(object), data)
{
#   ## Return the primary covariate
#   if (!(inherits(form, "formula"))) {
#     stop("\"Form\" must be a formula")
#   }
#   aux <- newgetCovariateFormula(form)
#   if (length(all.vars(aux)) > 0) {
#     eval(aux[[2]], object)
#   } else {
#     rep(1, dim(object)[1])
#   }
}

newgetGroups.data.frame <-
  ## Return the groups associated with object according to form for level
  function(object, form = formula(object), level, data, sep = "/")
{
  if (!missing(data)) {
    stop( "data argument to data.frame method for getGroups doesnt make sense" )
  }
  if (inherits(form, "formula")) {
    grpForm <- newgetGroupsFormula(form, asList = TRUE, sep = sep)
    if (is.null(grpForm)) {
      ## will use right hand side of form as the group formula
      grpForm <- newsplitFormula(asOneSidedFormula(form[[length(form)]]),
                              sep = sep)
      names(grpForm) <-
        unlist( lapply( grpForm, function(el) deparse( el[[ length(el) ]] ) ) )
    }
    if (any(unlist(lapply(grpForm,
#                          function(el) length(el[[length(el)]]))) != 1)) {
                          function(el) length(all.vars(el)))) != 1)) {
      stop("Invalid formula for groups")
    }
    form <- grpForm
  } else if (data.class(form) == "list") {
    if (!all(unlist(lapply(form, function(el) inherits(el, "formula"))))) {
      stop("Form must have all components as formulas")
    }
  } else {
    stop("Form can only be a formula, or a list of formulas")
  }
  vlist <- lapply(form,
                  function(x, dat, N) {
                    val <- eval(x[[length(x)]], dat)
                    if (length(val) == 1) {             # repeat groups
                      return(as.factor(rep(val, N)))
                    } else {
                      return(as.factor(val)[drop = TRUE])
                    }
                  }, dat = object, N = nrow(object))
  if (length(vlist) == 1) return(vlist[[1]]) # ignore level - only one choice
  ## make the list into a data frame with appropriate names
  value <- do.call("data.frame", vlist)
  if (missing(level)) return(value)
  if (is.character(level)) {
    nlevel <- match(level, names(vlist))
    if (any(aux <- is.na(nlevel))) {
      stop(paste("Level of", level[aux],"does not match formula \"",
		 deparse(form), "\""))
    }
  } else {
    nlevel <- as.numeric(level)
    if (any(aux <- is.na(match(nlevel, 1:ncol(value))))) {
      stop(paste("level of ", level[aux]," does not match formula \"",
	       deparse(form), "\""))
    }
  }
  if (length(nlevel) > 1)  return(value[, nlevel]) # multicolumn selection
  if (nlevel == 1)         return(value[, 1])     # no need to do more work
  value <- value[, 1:nlevel]
  val <- as.factor(do.call("paste", c(lapply(as.list(value),
					     as.character), sep = sep)))
  if (inherits(value[, 1], "ordered")) {
    value <- value[do.call("order", value),]
    aux <- unique(do.call("paste", c(lapply(as.list(value),
					    as.character), sep = sep)))
    return(ordered(val, aux))
  } else {
    return(ordered(val, unique(as.character(val))))
  }
}

getResponse.data.frame <-
  function(object, form = formula(object))
{
#   ## Return the response, the evaluation of the left hand side of a formula
#   ## on object
#   if (!(inherits(form, "formula") && (length(form) == 3))) {
#     stop("\"Form\" must be a two sided formula")
#   }
#   eval(form[[2]], object)
}

newgetGroupsFormula.default <-
  ## Return the formula(s) for the groups associated with object.
  ## The result is a one-sided formula unless asList is TRUE in which case
  ## it is a list of formulas, one for each level.
  function(object, asList = FALSE, sep = "/")
{
  form <- formula(object)
  if (!inherits(form, "formula")){
    stop("\"Form\" argument must be a formula")
  }
  form <- form[[length(form)]]
  if (!((length(form) == 3) && (form[[1]] == as.name("|")))) {
    ## no conditioning expression
    return(NULL)
  }
  ## val <- list( asOneSidedFormula( form[[ 3 ]] ) )
  val <- newsplitFormula(asOneSidedFormula(form[[3]]), sep = sep)
  names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
#  if (!missing(level)) {
#    if (length(level) == 1) {
#      return(val[[level]])
#    } else {
#      val <- val[level]
#    }
#  }
  if (asList) as.list(val)
  else as.formula(eval(parse(text = paste("~",  paste(names(val),
                               collapse = sep)))))
}

Names.formula <-
  function(object, data = list(), exclude = c("pi", "."), ...)
{
  if (!is.list(data)) { return(NULL) }  # no data to evaluate variable names
  allV <- all.vars(object)
  allV <- allV[is.na(match(allV, exclude))]

  if (length(allV) == 0) {
    if (attr(terms(object), "intercept")) { return("(Intercept)") }
    return(NULL)
  }

  if (any(is.na(match(allV, names(data))))) { return(NULL) }
  dimnames(model.matrix(object, model.frame(object, data)))[[2]]
}

Names.listForm <-
  function(object, data = list(), exclude = c("pi", "."), ...)
{
  pnames <- as.character(unlist(lapply(object, "[[", 2)))
  nams <- lapply(object, function(el, data, exclude) {
    Names(newgetCovariateFormula(el), data, exclude)
    }, data = data, exclude = exclude)
  if (is.null(nams[[1]])) return(NULL)
  val <- c()
  for(i in seq_along(object)) {
    if ((length(nams[[i]]) == 1) && (nams[[i]] == "(Intercept)")) {
      val <- c(val, pnames[i])
    } else {
      val <- c(val, paste(pnames[i], nams[[i]], sep = "."))
    }
  }
  val
}

newneedUpdate.default <-
  function(object)
{
#   val <- attr(object, "needUpdate")
#   if (is.null(val) || !val) FALSE
#   else TRUE
}

print.correlation <-
  ## Print only the lower triangle of a correlation matrix
  function(x, title = " Correlation:", rdig = 3, ...)
{
  p <- dim(x)[2]
  if (p > 1) {
    cat(title, "\n")
    ll <- lower.tri(x)
    x[ll] <- format(round(x[ll], digits = rdig))
    x[!ll] <- ""
    if (!is.null(colnames(x))) {
      colnames(x) <- abbreviate(colnames(x), minlength = rdig + 3)
    }
    print(x[-1,  - p, drop = FALSE], ..., quote = FALSE)
  }
  invisible(x)
}

###              Classes of positive-definite matrices

newpdConstruct <-
  ## a virtual constructor for these objects
  function(object, value, form, nam, data, ...) UseMethod("newpdConstruct")

pdFactor <-
  function(object) UseMethod("pdFactor")

pdMatrix <-
  ## extractor for the pd, correlation, or square-root factor matrix
  function(object, factor = FALSE) UseMethod("pdMatrix")

##*## pdMat - a virtual class of positive definite matrices

###*#  constructor for the virtual class

pdMat <-
  function(value = numeric(0), form = NULL, nam = NULL,
	   data = sys.frame(sys.parent()), pdClass = "pdSymm")
{
  if (inherits(value, "pdMat")) {	# nothing to construct
    pdClass <- class(value)
  }
  object <- numeric(0)
  class(object) <- unique(c(pdClass, "pdMat"))
  newpdConstruct(object, value, form, nam, data)
}

###*# Methods for local generics

corMatrix.pdMat <-
  function(object, ...)
{
  if (!newisInitialized(object)) {
    stop("Cannot access the matrix of uninitialized objects")
  }
  Var <- pdMatrix(object)
  if (length(unlist(dimnames(Var))) == 0) {
    aux <- paste("V", 1:(newDim(Var)[2]), sep = "")
    dimnames(Var) <- list(aux, aux)
  }
  dd <- dim(Var)
  dn <- dimnames(Var)
  stdDev <- sqrt(diag(Var))
  names(stdDev) <- colnames(Var)
  value <- array(t(Var/stdDev)/stdDev, dd, dn)
  attr(value, "stdDev") <- stdDev
  value
}

newpdConstruct.pdMat <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.frame(sys.parent()), ...)
{
  if (inherits(value, "pdMat")) {	# constructing from another pdMat
    if (length(form) == 0) {
      form <- formula(value)
    }
    if (length(nam) == 0) {
      nam <- Names(value)
    }
    if (newisInitialized(value)) {
      return(newpdConstruct(object, as.matrix(value), form, nam, data))
    } else {
      return(newpdConstruct(object, form = form, nam = nam, data = data))
    }
  }
  if (length(value) > 0) {
    if (inherits(value, "formula") || data.class(value) == "call") {
      ## constructing from a formula
      if (!is.null(form)) {
	warning("Ignoring argument \"form\"")
      }
      form <- formula(value)
      if (length(form) == 3) {          #two-sided case - nlme
        form <- list(form)
      }
    } else if (is.character(value)) {	# constructing from character array
      if (length(nam) > 0) {
	warning("Ignoring argument \"nam\"")
      }
      nam <- value
    } else if (is.matrix(value)) {	# constructing from a pd matrix
      vdim <- dim(value)
      if (length(vdim) != 2 || diff(vdim) != 0) {
        stop("\"value\" must be a square matrix")
      }
      if (length(unlist(vnam <- dimnames(value))) > 0) {
        vnam <- unique(unlist(vnam))
        if (length(vnam) != vdim[1]) {
          stop("dimnames of value must match or be NULL")
        }
        dimnames(value) <- list(vnam, vnam)
        if (length(nam) > 0) {          # check consistency
	  if (any(is.na(match(nam, vnam))) || any(is.na(match(vnam, nam)))) {
	    stop(paste("Names of \"value\" are not consistent",
		       "with \"nam\" argument"))
	  }
	  value <- value[nam, nam, drop = FALSE]
	} else {
	  nam <- vnam
	}
      }
      form <- form                      # avoid problems with lazy evaluation
      nam <- nam
      object <- chol((value + t(value))/2) # ensure it is positive-definite
      attr(object, "dimnames") <- NULL
      attr(object, "rank") <- NULL
    } else if (is.numeric(value)) {	# constructing from the parameter
      value <- as.numeric(value)
      attributes(value) <- attributes(object)
      object <- value
    } else if (data.class(value) == "list") {
      ## constructing from a list of two-sided formulae - nlme case
      if (!is.null(form)) {
	warning("Ignoring argument \"form\"")
      }
      form <- value
    } else {
      stop(paste(deparse(object), "is not a valid object for \"pdMat\""))
    }
  }

  if (!is.null(form)) {
    if (inherits(form, "formula") && length(form) == 3) {#two-sided case - nlme
      form <- list(form)
    }
    if (is.list(form)) {   # list of formulae
      if (any(!unlist(lapply(form,
                             function(el) {
                               inherits(el, "formula") && length(el) == 3
                             })))) {
        stop("All elements of \"form\" list must be two-sided formulas")
      }
      val <- list()
      for(i in seq_along(form)) {
        if (is.name(form[[i]][[2]])) {
          val <- c(val, list(form[[i]]))
        } else {
          val <- c(val, eval(parse(text = paste("list(",
            paste(paste(all.vars(form[[i]][[2]]), deparse(form[[i]][[3]]),
                        sep = "~"), collapse=","),")"))))
        }
      }
      form <- val
      class(form) <- "listForm"
      namesForm <- Names(form, data)
    } else {
      if (inherits(form, "formula")) {
        namesForm <- Names(asOneSidedFormula(form), data)
##        namesForm1 <- NULL
      } else {
        stop("\"form\" can only be a formula or a list of formulae")
      }
    }
    if (length(namesForm) > 0) {
      if (length(nam) == 0) {             # getting names from formula
        nam <- namesForm
      } else {				# checking consistency with names
        if (any(noMatch <- is.na(match(nam, namesForm)))) {
          err <- TRUE
          namCopy <- nam
          indNoMatch <- (1:length(nam))[noMatch]
          if (any(wch1 <- (nchar(nam, "c") > 12))) {
            ## possibly names with .(Intercept) in value
            wch1 <- substring(nam, nchar(nam, "c")-10) == "(Intercept)"
            if (any(wch1)) {
              namCopy[indNoMatch[wch1]] <-
                substring(nam[wch1], 1, nchar(nam[wch1], "c") - 12)
              noMatch[wch1] <- FALSE
              indNoMatch <- indNoMatch[!wch1]  # possibly not matched
            }
          }
          if (sum(noMatch) > 0) {
            ## still no matches - try adding .(Intercept)
            namCopy[indNoMatch] <-
              paste(namCopy[indNoMatch], "(Intercept)", sep = ".")
          }
          ## try matching modified value
          if (!any(is.na(match(namCopy, namesForm)))) {
            err <- FALSE
          }
          if (err) stop("\"form\" not consistent with \"nam\"")
        }
      }
    }
  }

  if (is.matrix(object)) {	# initialized as matrix, check consistency
    if (length(nam) > 0 && (length(nam) != dim(object)[2])) {
      stop(paste("Length of nam not consistent with dimensions",
		 "of initial value"))
    }
  }
  attr(object, "formula") <- form
  attr(object, "Dimnames") <- list(nam, nam)
  object
}

pdFactor.pdMat <-
  function(object)
{
  c(qr.R(qr(pdMatrix(object))))
}

pdMatrix.pdMat <-
  function(object, factor = FALSE)
{
  if (!newisInitialized(object)) {
    stop("Cannot access the matrix of uninitialized objects")
  }
  if (factor) {
    stop(paste("No default method for extracting the square",
               "root of a pdMat object"))
  } else {
    crossprod(pdMatrix(object, factor = TRUE))
  }
}

###*# Methods for standard generics

as.matrix.pdMat <-
  function(x, ...) pdMatrix(x)

newcoef.pdMat <-
  function(object, unconstrained = TRUE, ...)
{
  if (unconstrained || !newisInitialized(object)) {
    as.vector(object)
  } else {
    stop("Do not know how to obtain constrained coefficients")
  }
}

"coef<-.pdMat" <-
  function(object, ..., value)
{
  value <- as.numeric(value)
  if (newisInitialized(object)) {
    if (length(value) != length(object)) {
      stop("Cannot change the length of the parameter after initialization")
    }
  } else {
    return(newpdConstruct(object, value))
  }
  class(value) <- class(object)
  attributes(value) <- attributes(object)
  value
}

newDim.pdMat <-
  function(object, ...)
{
  if ((val <- length(Names(object))) > 0) {
    return(c(val, val))
  } else if (newisInitialized(object)) {
    return(dim(as.matrix(object)))
  }
  stop(paste("Cannot access the number of columns of",
	     "uninitialized objects without names."))
}

formula.pdMat <-
  function(x, asList, ...) eval(attr(x, "formula"))

#isInitialized.pdMat <- function(object) { newisInitialized.pdMat(object) }

newisInitialized.pdMat <-
  function(object)
{
  length(object) > 0
}

logDet.pdMat <-
  function(object, ...)
{
  if (!newisInitialized(object)) {
    stop(paste("Cannot extract the log of the determinant",
	       "from an uninitialized object"))
  }
  sum(log(svd(pdMatrix(object, factor = TRUE))$d))
}

"matrix<-.pdMat" <-
  function(object, value)
{
  value <- as.matrix(value)
  ## check for consistency of dimensions when object is initialized
  if (newisInitialized(object) && any(dim(value) != newDim(object))) {
    stop("Cannot change dimensions on an initialized pdMat object")
  }
  newpdConstruct(object, value)
}

Names.pdMat <-
  function(object, ...)
{
  return(as.character(attr(object, "Dimnames")[[2]])) 
}

"Names<-.pdMat" <-
  function(object, ..., value)
{
  if (is.null(value)) {
    attr(object, "Dimnames") <- NULL
    return(object)
  } else {
    value <- as.character(value)
    if (length(dn <- Names(object)) == 0) {
      if (newisInitialized(object)) {	# object is initialized without names
	if (length(value) != (aux <- newDim(object)[2])) {
	  stop(paste("Length of names should be", aux))
	}
      }
      attr(object, "Dimnames") <- list(value, value)
      return(object)
    }
    if (length(dn) != length(value)) {
      stop(paste("Length of names should be", length(dn)))
    }
    err <- FALSE
    if (any(noMatch <- is.na(match(value, dn)))) {
      err <- TRUE
      ## checking nlme case
      valueCopy <- value
      indNoMatch <- (1:length(value))[noMatch]
      nam1 <- value[noMatch]            # no matching names
      if (any(wch1 <- (nchar(nam1, "c") > 12))) {
        ## possibly names with .(Intercept) in value
        wch1 <- substring(nam1, nchar(nam1, "c")-10) == "(Intercept)"
        if (any(wch1)) {
          valueCopy[indNoMatch[wch1]] <-
            substring(nam1[wch1], 1, nchar(nam1[wch1], "c") - 12)
          noMatch[wch1] <- FALSE
          indNoMatch <- indNoMatch[!wch1]  # possibly not matched
        }
      }
      if (sum(noMatch) > 0) {
        ## still no matches - try adding .(Intercept)
        valueCopy[indNoMatch] <-
          paste(valueCopy[indNoMatch], "(Intercept)", sep = ".")
      }
      ## try matching modified value
      indMatch <- match(valueCopy, dn)
      if (!any(is.na(indMatch))) {      # all match
        attr(object, "Dimnames") <- list(value, value)
        if ((length(indMatch)) > 1 && any(diff(indMatch) != 1) &&
            newisInitialized(object)) { # permutation
          auxMat <- as.matrix(object)[indMatch, indMatch, drop = FALSE]
          dimnames(auxMat) <- list(value, value)
          return(newpdConstruct(object, auxMat))
        }
        return(object)
      }
    }
    if (err) {
      stop(paste("Names being assigned do not correspond to a permutation",
                 "of previous names", sep = "\n"))
    }
    indMatch <- match(value, dn)
    if ((length(indMatch) == 1) || all(diff(indMatch) == 1)) {
      return(object)
    }
    ## must be a permutation of names
    attr(object, "Dimnames") <- list(value, value)
    if (newisInitialized(object)) {
      auxMat <- as.matrix(object)[indMatch, indMatch, drop = FALSE]
      dimnames(auxMat) <- list(value, value)
      return(newpdConstruct(object, auxMat))
    }
    object
  }
}

plot.pdMat <-
  function(x, nseg = 50, levels = 1, center = rep(0, length(stdDev)),
	   additional, ...)
{
  corr <- corMatrix(x)
  stdDev <- attr(corr, "stdDev")
  attr(corr, "stdDev") <- NULL
  assign(".corr", corr)
  assign(".angles", seq(-pi, pi, length = nseg + 1))
  assign(".cosines", cos(.angles))
  nlev <- length(levels)
  dataMat <- array(aperm(outer(rbind(-stdDev, stdDev), levels), c(1, 3, 2)),
		   dim = c(nlev * 2, length(stdDev)),
		   dimnames = list(NULL, names(stdDev)))
  groups <- rep(1:nlev, rep(2, nlev))
  dataMat <- t(t(dataMat) + center)
  if (!missing(additional)) {
    additional <- as.matrix(additional)
    dataMat <- rbind(dataMat, additional)
    groups <- c(groups, rep(0, nrow(additional)))
  }
  splom(~ dataMat, panel = function(x, y, subscripts, groups, ...) {
    groups <- groups[subscripts]	# should be a no-op but
    if (any(g0 <- groups == 0)) {	# plot as points
      panel.xyplot(x[g0], y[g0], ..., type = "p")
    }
    g1 <- groups == 1			# plot the center points
    panel.xyplot(mean(x[g1]), mean(y[g1]), ..., type = "p", pch = 3)
    p <- ncol(.corr)
    laggedCos <- cos(.angles + acos(.corr[round(mean(x[g1])*p + 0.5),
					  round(mean(y[g1])*p + 0.5)]))
    xylist <- lapply(split(data.frame(x = x[!g0], y = y[!g0]), groups[!g0]),
		     function(el, lagged) {
		       if (nrow(el) != 2) {
			 stop("x-y data to splom got botched somehow")
		       }
		       sumDif <- array(c(1,1,1,-1)/2, c(2,2)) %*% as.matrix(el)
		       list(x = sumDif[1,1] + .cosines * sumDif[2,1],
			    y = sumDif[1,2] + lagged * sumDif[2,2])
		     }, lagged = laggedCos)
    gg <- rep(seq_along(xylist), rep(length(.angles), length(xylist)))
    panel.superpose(unlist(lapply(xylist, "[[", "x")),
		    unlist(lapply(xylist, "[[", "y")),
		    subscripts = seq_along(gg), groups = gg, ..., type = "l")
  }, subscripts = TRUE, groups = groups)
}

print.pdMat <-
  function(x, ...)
{
  if (newisInitialized(x)) {
    cat("Positive definite matrix structure of class", class(x)[1], "representing\n")
    print(as.matrix(x), ...)
  } else {
    cat("Uninitialized positive definite matrix structure of class ", class(x)[1],
	".\n", sep = "")
  }
  invisible(x)
}

print.summary.pdMat <-
  function(x, sigma = 1, rdig = 3, Level = NULL, resid = FALSE, ...)
  ## resid = TRUE causes an extra row to be added
{
  if (!is.list(x)) {
    if (!(is.null(form <- attr(x, "formula")))) {
      cat(paste(" Formula: "))
      if (inherits(form, "formula")) {
        cat(deparse(form))
        if (!is.null(Level)) { cat( paste( " |", Level ) ) }
      } else {
        if (length(form) == 1) {
          cat(deparse(form[[1]]))
          if (!is.null(Level)) { cat( paste( " |", Level ) ) }
        } else {
          cat(deparse(lapply(form,
                             function(el) as.name(deparse(el)))))
          cat("\n Level:", Level)
        }
      }
      cat( "\n" )
    }
    if (ncol(x) == 1) {
      if (resid) {
        print(array(sigma * c(attr(x, "stdDev"), 1), c(1, 2),
                    list("StdDev:",
                         c(names(attr(x, "stdDev")), "Residual"))), ... )
      } else {
        print(array(sigma * attr(x, "stdDev"), c(1,1),
                    list("StdDev:", names(attr(x, "stdDev")))), ... )
      }
    } else {
      cat(paste(" Structure: ", attr(x, "structName"), "\n", sep = ""))
      if (attr(x, "noCorrelation") | (1 >= (p <- dim(x)[2]))) {
        if (resid) {
          print(array(sigma * c(attr(x, "stdDev"), 1), c(1, p + 1),
                      list("StdDev:",
                           c(names(attr(x, "stdDev")), "Residual"))), ...)
        } else {
          print(array(sigma * attr(x, "stdDev"), c(1, p),
                      list("StdDev:", names(attr(x, "stdDev")))), ...)
        }
      } else {                          # we essentially do print.correlation here
        ll <- lower.tri(x)
        stdDev <- attr(x, "stdDev")
        x[ll] <- format(round(x[ll], digits = rdig), ...)
        x[!ll] <- ""
        xx <- array("", dim(x),
                    list(names(attr(x, "stdDev")),
                         c("StdDev", "Corr", rep("", p - 2))))
        xx[, 1] <- format(sigma * attr(x, "stdDev"))
        xx[-1, -1] <- x[ -1, -p ]
        if (!is.null(colnames(x))) {
          xx[1, -1] <- abbreviate(colnames(x)[ -p ], minlength = rdig + 3)
        }
        if (resid) {
          x <- array("", dim(xx) + c(1, 0),
                     list(c(rownames(xx), "Residual"), colnames(xx)))
          x[ 1:p, ] <- xx
          x[ , 1 ] <- format(sigma * c(stdDev, 1))
          xx <- x
        }
        print( xx, ..., quote = FALSE )
      }
    }
  } else {				# composite structure
    cat(paste(" Composite Structure: ", attr(x, "structName"), "\n", sep =""))
    elName <- attr(x, "elementName")
    compNames <- names(x)
    for (i in seq_along(x)) {
      cat(paste("\n ", elName, " ", i, ": ", compNames[i], "\n", sep = ""))
      print.summary.pdMat(x[[i]], sigma = sigma, Level = Level,
                          resid = resid && (i == length(x)), ...)
    }
  }
  invisible(x)
}

solve.pdMat <-
  function(a, b, ...)
{
  if (!newisInitialized(a)) {
    stop("Cannot get the inverse of an uninitialized object")
  }
  matrix(a) <- solve(as.matrix(a))
  a
}

summary.pdMat <-
  function(object, structName = class(object)[1], noCorrelation = FALSE, ...)
{
  #cat('MY OUTPUT - summary.pdMat is running from lmeTogether\n')
  if (newisInitialized(object)) {
    value <- corMatrix(object)
    attr(value, "structName") <- structName
    attr(value, "noCorrelation") <- noCorrelation
    attr(value, "formula") <- formula(object)
    class(value) <- "summary.pdMat"
    value
  } else {
    object
  }
}

"[.pdMat" <-
  function(x, i, j, drop = TRUE)
{
#   xx <- x
#   x <- as.matrix(x)
#   if (missing(i)) li <- 0
#   else li <- length(i)
#   if (missing(j)) lj <- 0
#   else lj <- length(j)
# 
#   if ((li + lj == 0) ||
#       (li == lj) && ((mode(i) == mode(j)) && all(i == j))) {
#     drop <- FALSE			# even for a 1 by 1 submatrix,
# 					# you want it to be a matrix
#     newpdConstruct(xx, NextMethod())
#   } else {
#     NextMethod()
#   }
}

"[<-.pdMat" <-
  function(x, i, j, value)
{
#   xx <- x
#   x <- as.matrix(x)
#   newpdConstruct(xx, NextMethod())
}

##*## Classes that substitute for (i.e. inherit from) pdMat

###*# pdSymm - a class of general pd matrices

####* Constructor

pdSymm <-
  ## Constructor for the pdSymm class
  function(value = numeric(0), form = NULL, nam = NULL, data = parent.frame())
{
#   object <- numeric(0)
#   class(object) <- c("pdSymm", "pdMat")
#   newpdConstruct(object, value, form, nam, data)
}

####* Methods for local generics

newpdConstruct.pdSymm <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.frame(sys.parent()), ...)
{
#   val <- NextMethod()
#   if (length(val) == 0) {               # uninitialized object
#     class(val) <- c("pdSymm", "pdMat")
#     return(val)
#   }
# 
#   if (is.matrix(val)) {
#     vald <- svd(val, nu = 0)
#     object <- vald$v %*% (log(vald$d) * t(vald$v))
#     value <- object[row(object) <= col(object)]
#     attributes(value) <- attributes(val)[names(attributes(val)) !=  "dim"]
#     class(value) <- c("pdSymm", "pdMat")
#     return(value)
#   }
#   Ncol <- round((sqrt(8*length(val) + 1) - 1)/2)
#   if (length(val) != round((Ncol * (Ncol + 1))/2)) {
#     stop(paste("An object of length", length(val),
# 	       "does not match the required parameter size"))
#   }
#   class(val) <- c("pdSymm", "pdMat")
#   val
}

pdFactor.pdSymm <-
  function(object)
{
  Ncol <- round((-1 + sqrt(1 + 8 * length(object))) / 2)
  .C("matrixLog_pd",
     Factor = double(Ncol * Ncol),
     as.integer(Ncol),
     as.double(object))$Factor
}

pdMatrix.pdSymm <-
  function(object, factor = FALSE)
{
  if (!newisInitialized(object)) {
    stop("Cannot extract matrix from an uninitialized object")
  }
  if (factor) {
    Ncol <- newDim(object)[2]
    value <- array(pdFactor(object), c(Ncol, Ncol), attr(object, "Dimnames"))
    attr(value, "logDet") <- sum(log(abs(svd(value)$d)))
    value
  } else {
    NextMethod()
  }
}

####* Methods for standard generics

newcoef.pdSymm <-
  function(object, unconstrained = TRUE, ...)
{
  if (unconstrained || !newisInitialized(object)) NextMethod()
  else {				# upper triangular elements
    val <- as.matrix(object)
    aN <- Names(object)
    aN1 <- paste("cov(", aN, sep ="")
    aN2 <- paste(aN, ")", sep ="")
    aNmat <- t(outer(aN1, aN2, paste, sep = ","))
    aNmat[row(aNmat) == col(aNmat)] <- paste("var(",aN,")",sep="")
    val <- val[row(val) <= col(val)]
    names(val) <- aNmat[row(aNmat) <= col(aNmat)]
    val
  }
}

newDim.pdSymm <-
  function(object, ...)
{
  if (newisInitialized(object)) {
    val <- round((sqrt(8*length(object) + 1) - 1)/2)
    c(val, val)
  } else {
    NextMethod()
  }
}

logDet.pdSymm <-
  function(object, ...)
{
  if (!newisInitialized(object)) {
    stop(paste("Cannot extract the log of the determinant",
	       "from an uninitialized object"))
  }
  attr(pdMatrix(object, factor = TRUE), "logDet")
}

solve.pdSymm <-
  function(a, b, ...)
{
  if (!newisInitialized(a)) {
    stop("Cannot extract the inverse from an uninitialized object")
  }
  coef(a) <- -coef(a, TRUE)
  a
}

summary.pdSymm <-
  function(object,
	   structName = "General positive-definite", ...)
{
  cat('MY OUTPUT - summary.pdSymm is running from lmeTogether\n')
  summary.pdMat(object, structName)
}

### No need to implement other methods as the methods for pdMat
### are sufficient.

###*# pdLogChol - a general positive definite structure parameterized
###   by the non-zero elements of the Cholesky factor with the diagonal
###   elements given in the logarithm scale.

####* Constructor

pdLogChol <-
  ## Constructor for the pdLogChol class
  function(value = numeric(0), form = NULL, nam = NULL, data = sys.parent())
{
  object <- numeric(0)
  class(object) <- c("pdLogChol", "pdMat")
  newpdConstruct(object, value, form, nam, data)
}

####* Methods for local generics

newpdConstruct.pdLogChol <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.parent(), ...)
{
  val <- newpdConstruct.pdMat(object, value, form, nam, data)
  if (length(val) == 0) {               # uninitialized object
    class(val) <- c("pdLogChol", "pdSymm", "pdMat")
    return(val)
  }
  if (is.matrix(val)) {
    value <- c(log(diag(val)), val[row(val) < col(val)])
    attributes(value) <- attributes(val)[names(attributes(val)) != "dim"]
    class(value) <- c("pdLogChol", "pdSymm", "pdMat")
    return(value)
  }
  Ncol <- round((sqrt(8*length(val) + 1) - 1)/2)
  if (length(val) != round((Ncol * (Ncol + 1))/2)) {
    stop(paste("An object of length", length(val),
	       "does not match a Cholesky factor"))
  }
  class(val) <- c("pdLogChol", "pdSymm", "pdMat")
  val
}

pdFactor.pdLogChol <-
  function(object)
{
  Ncol <- round((-1 + sqrt(1 + 8 * length(object))) / 2)
  .C("logChol_pd",
     Factor = double(Ncol * Ncol),
     as.integer(Ncol),
     as.double(object))$Factor
}

####* Methods for standard generics

solve.pdLogChol <-
  function(a, b, ...)
{
  if (!newisInitialized(a)) {
    stop("Cannot get the inverse of an uninitialized object")
  }
  Ncol <- (-1 + sqrt(1 + 8 * length(a))) / 2
#  val <- array(.Fortran("dbksl",
# 			as.double(pdFactor(a)),
# 			as.integer(Ncol),
# 			as.integer(Ncol),
# 			val = as.double(diag(Ncol)),
# 			as.integer(Ncol),
# 			integer(1))[["val"]], c(Ncol, Ncol))
#  val <- qr(t(val))$qr
  val <- qr(t(solve(pdMatrix(a, factor = TRUE))))$qr
  val <- sign(diag(val)) * val
  coef(a) <- c(log(diag(val)), val[c(row(val) < col(val))])
  a
}

summary.pdLogChol <-
  function(object, structName =
           "General positive-definite, Log-Cholesky parametrization",
           ...)
{
  summary.pdMat(object, structName)
}


#### No need to implement other methods as the methods for pdMat
#### are sufficient.

#newpdConstruct.pdSymm <- newpdConstruct.pdMatrixLog    #default parametrization

####*# pdNatural - a general positive definite structure parameterized
####   by the log of the square root of the diagonal elements and the
####   generalized logit of the correlations. This is NOT an unrestricted
####   parametrization

####* Constructor

pdNatural <-
  ## Constructor for the pdNatural class
  function(value = numeric(0), form = NULL, nam = NULL, data = sys.frame(sys.parent()))
{
  object <- numeric(0)
  class(object) <- c("pdNatural", "pdMat")
  newpdConstruct(object, value, form, nam, data)
}

####* Methods for local generics

newpdConstruct.pdNatural <-
  function(object, value = numeric(0), form = formula(object),
	   nam = Names(object), data = sys.frame(sys.parent()), ...)
{
  val <- newpdConstruct.pdMat(object, value, form, nam, data)
  if (length(val) == 0) {               # uninitiliazed object
    class(val) <- c("pdNatural", "pdMat")
    return(val)
  }
  if (is.matrix(val)) {
    q <- ncol(val)
    if (q > 1) {
      aux <- crossprod(val)
      stdDev <- sqrt(diag(aux))
      aux <- t(aux/stdDev)/stdDev
      aux <- aux[row(aux) > col(aux)]
      value <- c(log(stdDev), log((aux + 1)/(1 - aux)))
    } else {
      value <- log(val)
    }
    attributes(value) <- attributes(val)[names(attributes(val)) != "dim"]
    class(value) <- c("pdNatural", "pdMat")
    return(value)
  }
  Ncol <- round((sqrt(8*length(val) + 1) - 1)/2)
  if (length(val) != round((Ncol * (Ncol + 1))/2)) {
    stop(paste("An object of length", length(val),
	       "does not match the required parameter size"))
  }
  class(val) <- c("pdNatural", "pdMat")
  val
}

pdFactor.pdNatural <-
  function(object)
{
  Ncol <- round((-1 + sqrt(1 + 8 * length(object))) / 2)
  .C("natural_pd",
     Factor = double(Ncol * Ncol),
     as.integer(Ncol),
     as.double(object))$Factor
}

pdMatrix.pdNatural <-
  function(object, factor = FALSE)
{
  if (!newisInitialized(object)) {
    stop("Cannot extract matrix from an uninitialized object")
  }
  if (factor) {
    Ncol <- newDim(object)[2]
    value <- array(pdFactor(object), c(Ncol, Ncol), attr(object, "Dimnames"))
    attr(value, "logDet") <- sum(log(diag(value)))
    value
  } else {
    NextMethod()
  }
}

####* Methods for standard generics

newcoef.pdNatural <-
  function(object, unconstrained = TRUE, ...)
{
  if (unconstrained || !newisInitialized(object)) NextMethod()
  else {				# standard deviations and correlations
    Ncol <- round((-1 + sqrt(1 + 8 * length(object))) / 2)
    val <- exp(as.vector(object))
    aux <- val[-(1:Ncol)]
    val[-(1:Ncol)] <- (aux - 1) / (aux + 1)
    aN <- Names(object)
    aNmat <- t(outer(aN, aN, paste, sep = ","))
    names(val) <- c(paste("sd(",aN,")", sep = ""),
		    if (Ncol > 1) {
		      paste("cor(", aNmat[row(aNmat) > col(aNmat)],")",sep="")
		    })
    val
  }
}

newDim.pdNatural <-
  function(object, ...)
{
  if (newisInitialized(object)) {
    val <- round((sqrt(8*length(object) + 1) - 1)/2)
    c(val, val)
  } else {
    NextMethod()
  }
}

logDet.pdNatural <-
  function(object, ...)
{
  if (!newisInitialized(object)) {
    stop(paste("Cannot extract the log of the determinant",
	       "from an uninitialized object"))
  }
  attr(pdMatrix(object, factor = TRUE), "logDet")
}


solve.pdNatural <-
  function(a, b, ...)
{
  if (!newisInitialized(a)) {
    stop("Cannot get the inverse of an uninitialized object")
  }
  Ncol <- round((-1 + sqrt(1 + 8 * length(a))) / 2)
  if (Ncol > 1) {
#     val <- array(.Fortran("dbksl",
# 			  as.double(pdFactor(a)),
# 			  as.integer(Ncol),
# 			  as.integer(Ncol),
# 			  val = as.double(diag(Ncol)),
# 			  as.integer(Ncol),
# 			  integer(1))[["val"]], c(Ncol, Ncol))
    val <- solve(pdMatrix(a, factor = TRUE))
    val <- val %*% t(val)
    stdDev <- sqrt(diag(val))
    val <- t(val/stdDev)/stdDev
    val <- val[row(val) > col(val)]
    coef(a) <- c(log(stdDev), log((val + 1)/(1 - val)))
  } else {
    coef(a) <- -coef(a)
  }
  a
}

summary.pdNatural <-
  function(object,
	   structName = "General positive-definite, Natural parametrization",
           ...)
{
  summary.pdMat(object, structName)
}

###      Methods for the class of random-effects structures.

##*## Generics that should be implemented for any reStruct class

###*# Constructor

reStruct <-
  function(object, pdClass = "pdLogChol", REML = FALSE, data = sys.frame(sys.parent()))
{
  ## object can be:
  ## 1) a named list of formulas or pdMats with grouping factors as names
  ##    (assume same order of nesting as order of names)
  ## 2) a formula of the form ~ x | g or ~ x | g1/g2/../gn
  ## 3) a list of formulas like ~x | g
  ## 4) a formula like ~x, a pdMat object, or a list of such
  ##    formulas or objects . In this case, the data used to
  ##    initialize the reStruct will be required to inherit from class
  ##    "groupedData"
  ## 5) another reStruct object
  ## parametrization specifies the pdMat constructor to be used for all
  ## formulas used in object

  if (inherits(object, "reStruct")) {	# little to do, return object
    if (!missing(REML)) attr(object, "settings")[1] <- as.integer(REML)
    object[] <- lapply(object,
		       function(el, data) {
			 pdMat(el, data = data)
		       }, data = data)
    return(object)
  }
  plen <- NULL
  if (inherits(object, "formula")) {	# given as a formula
    if (is.null(grpForm <- newgetGroupsFormula(object, asList = TRUE))) {
      cat('MY OUTPUT - object is null from reStruct is running from lmeTogether\n')
      object <- list( object )
    } else {
      if (length(object) == 3) {        # nlme type of formula
        object <-
          eval(parse(text = paste(deparse(newgetResponseFormula(object)[[2]]),
                       deparse(newgetCovariateFormula(object)[[2]], width.cutoff=500),
                     sep = "~")))
      } else {
        object <- newgetCovariateFormula(object)
      }
      object <- rep( list(object), length( grpForm ) )
      names( object ) <- names( grpForm )
    }
  } else if (inherits(object, "pdMat")) { # single group, as pdMat
    if (is.null(formula(object))) {
      stop("pdMat element must have a formula")
    }
    object <- list(object)
  } else {
    if (data.class(object) != "list") {
      stop("Object must be a list or a formula")
    }
    ## checking if nlme-type list - unnamed list of 2-sided formulas
    if (is.null(names(object)) &&
        all(unlist(lapply(object, function(el) {
          inherits(el, "formula") && length(el) == 3})))) {
      object <- list(object)
    } else {
      ## checking if elements are valid
      object <- lapply(object,
                       function(el) {
                         if (inherits(el, "pdMat")) {
                           if (is.null(formula(el))) {
                             stop("pdMat elements must have a formula")
                           }
                           return(el)
                         }
                         if (inherits(el, "formula")) {
                           grpForm <- newgetGroupsFormula(el)
                           if (!is.null(grpForm)) {
                             el <- newgetCovariateFormula(el)
                             attr(el, "grpName") <- deparse(grpForm[[2]])
                           }
                           return(el)
                         } else {
                           if (data.class(el) == "list" &&
                               all(unlist(lapply(el, function(el1) {
                                 inherits(el1, "formula") && length(el1) == 3
                               })))) { return(el) }
                           else {
                 stop("Elements in object must be formulas or pdMat objects")
                           }
                         }
		     })
    }
    if (is.null(namObj <- names(object))) {
      namObj <- rep("", length(object))
    }
    aux <- unlist(lapply(object,
			 function(el) {
			   if (inherits(el, "formula") &&
			       !is.null(attr(el, "grpName"))) {
			     attr(el, "grpName")
			   } else ""
			 }))
    auxNam <- namObj == ""
    if (any(auxNam)) {
      namObj[auxNam] <- aux[auxNam]
    }
    names(object) <- namObj
  }

  ## converting elements in object to pdMat objects
  object <- lapply(object,
		   function(el, pdClass, data) {
#                     if (data.class(el) == "pdSymm")
#                       warning("class pdSymm may cause problems if using analytic gradients")
		     pdMat(el, pdClass = pdClass, data = data)
		   }, pdClass = pdClass, data = data)

  object <- rev(object)			# inner to outer groups
  if (all(unlist(lapply(object, newisInitialized)))) {
    plen <- unlist(lapply(object, function(el) length(coef(el))))
  }
  pC <- unlist(lapply(object, data.class))
  pC <- match(pC, c("pdSymm", "pdDiag", "pdIdent", "pdCompSymm",
                    "pdLogChol"), 0) - 1
#  if (any(pC == -1)) {                 # multiple nesting
#    pC <- -1
#  }
  ## at this point, always require asDelta = TRUE and gradHess = 0
  attr(object, "settings") <- c(as.integer(REML), 1, 0, pC)
  attr(object, "plen") <- plen
  class(object) <- "reStruct"
  object
}

###*# Methods for pdMat generics

corMatrix.reStruct <-
  function(object, ...)
{
  if (!newisInitialized(object)) {
    stop("Cannot access the matrix of uninitialized objects")
  }
  as.list(rev(lapply(object, corMatrix)))
}

pdFactor.reStruct <-
  function(object)
{
  unlist(lapply(object, pdFactor))
}

pdMatrix.reStruct <-
  function(object, factor = FALSE)
{
  if (!newisInitialized(object)) {
    stop("Cannot access the matrix of uninitialized objects")
  }
  as.list(rev(lapply(object, pdMatrix, factor)))
}

###*# Methods for standard generics

as.matrix.reStruct <-
  function(x, ...) pdMatrix(x)

newcoef.reStruct <-
  function(object, unconstrained = TRUE, ...)
{
  unlist(lapply(object, coef, unconstrained))
}

# "coef<-.reStruct" <-
#   function(object, ..., value)
# {
#   if (is.null(plen <- attr(object, "plen"))) {
#     stop(paste("Cannot change the parameter when",
# 	       "length of parameters is undefined"))
#   }
#   if (length(value) != sum(plen)) {
#     stop("Cannot change parameter length of initialized objects")
#   }
#   ends <- cumsum(plen)
#   starts <- 1 + c(0, ends[-length(ends)])
#   #print(object)
#   for (i in seq_along(object)) {
#     coef(object[[i]]) <- value[(starts[i]):(ends[i])]
#   }
#   object
# }

"coef<-.reStruct" <-
  function(object, ..., value)
{
  if (is.null(plen <- attr(object, "plen"))) {
    stop(paste("Cannot change the parameter when",
	       "length of parameters is undefined"))
  }
  if (length(value) != sum(plen)) {
    stop("Cannot change parameter length of initialized objects")
  }
  ends <- cumsum(plen)
  starts <- 1 + c(0, ends[-length(ends)])
  for (i in seq_along(object)) {
    coef(object[[i]]) <- value[(starts[i]):(ends[i])]
  }
  object
}

newformula.reStruct <-
  function(x, asList = FALSE, ...)
{
  as.list(lapply(x, formula, asList))
}

newgetGroupsFormula.reStruct <-
  function(object, asList = FALSE, sep)
{
  if (is.null(val <- rev(formula(object)))) {
    stop("Can not extract groups formula without a formula")
  }
  if (is.null(nVal <- names(val))) return(NULL)
  if (asList) {
    for(i in nVal) {
      val[[i]] <- eval(parse(text = paste("~",i)))
    }
  } else {
    val <- eval(parse(text = paste("~",paste(nVal, collapse = "/"))))
  }
  val
}

newisInitialized.reStruct <-
  function(object)
{
 all(unlist(lapply(object, newisInitialized)))
}

newInitialize.reStruct <-
  function(object, data, conLin, control = list(niterEM = 20), ...)
{

  ## initialize reStruct object, possibly getting initial estimates
  seqO <- seq_along(object)
  ## check if names are defined
  lNams <- unlist(lapply(object, function(el) length(Names(el)))) == 0
  if (any(lNams)) {			# need to resolve formula names
    aux <- seqO[lNams]
    object[aux] <- lapply(object[aux],
			  function(el, data) {
			    newpdConstruct(el, el, data = data)
			  }, data = data)
  }

  ## obtaining the parameters mapping
  plen <- unlist(lapply(object, function(el)
			{

			  if (newisInitialized(el)) {
			    length(coef(el))			    
			  } else {
			    matrix(el) <- diag(length(Names(el)))
			    length(coef(el))
  			  }
			}))

  if (!all(plen > 0)) {
    stop("All elements of a reStruct object must have a non-zero size")
  }
  attr(object, "plen") <- plen

  ## checking initialization
  isIni <- unlist(lapply(object, newisInitialized))
  if (!all(isIni)) {			# needs initialization
    dims <- conLin$dims
    Q <- dims$Q
    qvec <- dims$qvec[1:Q]
    auxInit <-
      lapply(split(0.375^2 * apply((conLin$Xy[, 1:sum(qvec), drop = FALSE])^2,
	     2, sum)/ rep(dims$ngrps[1:Q], qvec), rep(1:Q, qvec)),
	     function(x) diag(x, length(x)))
  }
  for(i in seqO) {
    if (isIni[i]) {
      object[[i]] <- solve(object[[i]])	#working with precisions
    } else {
      matrix(object[[i]]) <- auxInit[[i]]
    }
    NULL
  }
  newMEEM(object, conLin, control$niterEM) # refine initial estimates with EM
}

newlogDet.reStruct <-
  function(object, ...)
{
  unlist(lapply(object, logDet))
}

logLik.reStruct <-
  function(object, conLin, ...)
{
  if(any(!is.finite(conLin$Xy))) return(-Inf)
  .C("mixed_loglik",
     as.double(conLin$Xy),
     as.integer(unlist(conLin$dims)),
     as.double(pdFactor(object)),
     as.integer(attr(object, "settings")),
     loglik = double(1),
     double(1))$loglik
}

"matrix<-.reStruct" <-
  function(object, value)
{
  if (data.class(value) != "list") value <- list(value)
  if (length(value) != length(object)) {
    stop("Cannot change the length of object")
  }
  value <- rev(value)                   # same order as object
  for(i in seq_along(object)) {
    matrix(object[[i]]) <- value[[i]]
  }
  object
}

# model.matrix.reStruct <-
#   function(object, data, contrast = NULL, ...) { model.matrix.reStruct(object,data,contrast, ...) }

model.matrix.reStruct <-
  function(object, data, contrast = NULL, ...)
{
  if (is.null(form <- formula(object, asList = TRUE))) {
    stop("Cannot extract model matrix without formula")
  }
  form1 <- newasOneFormula(form)
  if (length(form1) > 0) {
    data <- model.frame(form1, data = data)
  } else {
    data <- data.frame("(Intercept)" = rep(1, nrow(data)))
  }
  any2list <- function( object, data, contrast ) {
    form2list <- function(form, data, contrast) {
      if (length(newasOneFormula( form )) == 0) {# the ~ 1 case
        return(list("(Intercept)" = rep(1, dim(data)[1])))
      }
      as.data.frame(unclass(model.matrix(form,
                                         model.frame(form, data),
                                         contrast)))
    }
    if (inherits( object, "formula" )) {
      return( form2list( object, data, contrast ) )
    }
    if (is.list( object ) ) {
      return( unlist(lapply(object, form2list, data = data, contrast = contrast),
                     recursive = FALSE ) )
    }
    return( NULL)
  }
  value <- as.list(lapply(form, any2list,
                          data = data, contrast = contrast))
  ## save the contrasts currently in effect for later predictions
  contr <- as.list(lapply( as.data.frame(data), function(x)
                  if( inherits( x, "factor" ) &&
                     length(levels(x)) > 1) contrasts(x) else NULL ))
  contr[names(contrast)] <- contrast

  ncols <- as.vector(unlist(lapply(value, length)))
  nams <- if (length(value) == 1) {
    names(value[[1]])
  } else {
    paste(rep(names(value), ncols), unlist(lapply(value, names)), sep = ".")
  }
  val <- matrix(unlist(value), nrow = nrow(data),
                dimnames = list(row.names(data), nams))
  attr(val, "ncols") <- ncols
  attr(val, "nams") <- as.list(lapply(value, names))
  attr(val, "contr") <- contr
  val
}

Names.reStruct <-
    function(object, ...)
{
    as.list(lapply(object, Names))
}

"Names<-.reStruct" <-
  function(object, ..., value)
{
  if (length(object) != length(value)) {
    stop("Incompatible lengths for object names")
  }
  for(i in seq_along(object)) {
    Names(object[[i]]) <- value[[i]]
  }
  object
}

needUpdate.reStruct <-
  function(object) FALSE

print.reStruct <-
  function(x, sigma = 1, reEstimates, verbose = FALSE, ...)
{
  ox <- x
  if (newisInitialized(x)) {
    nobj <- length(x)
    if (is.null(namx <- names(x))) names(x) <- nobj:1
    aux <- t(array(rep(names(x), nobj), c(nobj, nobj)))
    aux[lower.tri(aux)] <- ""
    x[] <- rev(x)
    names(x) <-
      rev(apply(aux, 1, function(x) paste(x[x != ""], collapse = " %in% ")))
    cat("Random effects:\n")
    for(i in seq_along(x)) {
      print(summary(x[[i]]), sigma, Level = names(x)[i],
            resid = (i == length(x)), ...)
      if (verbose) {
	cat("Random effects estimates:\n")
	print(reEstimates[[i]])
      }
      cat("\n")
    }
  } else {
    cat("Uninitialized random effects structure\n")
  }
  invisible(ox)
}

newrecalc.reStruct <-
  function(object, conLin, ...)
{
  conLin[["logLik"]] <- conLin[["logLik"]] + logLik(object, conLin)
  conLin
}

solve.reStruct <-
  function(a, b, ...)
{
  a[] <- lapply(a, solve)
  a
}

summary.reStruct <- function(object, ...) object

update.reStruct <-
  function(object, data, ...)
{
  object
}

"[.reStruct" <-
  function(x, ...)
{
  val <- NextMethod()
  if (length(val)) class(val) <- "reStruct"
  val
}

###              Classes of variance functions

##*## Generics that should be implemented for any varFunc class

newvarWeights <-
  ## Calculates the weights of the variance function
  function(object) UseMethod("varWeights")

##*## varFunc - a virtual class of variance functions

###*# Constructor

newvarFunc <-
  ## Can take as argument either a varFunc object, in which case it does
  ## nothing, a formula or a character string , in which case it
  ## calls varFixed
  function(object)
{
  if(is.null(object)) return(object)	# NULL object - no varFunc structure
  if (inherits(object, "varFunc")) {
    ## constructing from another varFunc object
    return(object)
  }
  if (inherits(object, "formula") || is.character(object)) {
    ## constructing from a formula of the form ~ x
    return(varFixed(asOneSidedFormula(object)))
  }

  stop(paste("Can only construct varFunc object from another varFunc",
	     "object, a formula, or a character string"))
}


###*# Methods for local generics

varWeights.varFunc <-
  function(object) attr(object, "weights")


###   Miscellaneous methods that must be defined last in the library

# newBIC.lme <- newBIC.gls <- newBIC.lm
# 
# comparePred.lme <- comparePred.gls
# 
# getData.nlme <- getData.gnls
# 
# getData.lme <- getData.gls <- getData.nls
# 
# qqnorm.gls <- qqnorm.lm <- qqnorm.nls
# 
# plot.lme <- plot.nls
# 
# fitted.gnls <- fitted.gls
# 
# residuals.gnls <- residuals.gls

vcov.gls <- function (object, ...) object$varBeta

vcov.lme <- function (object, ...) object$varFix

.onUnload <- function(libpath)
    library.dynam.unload("nlme", libpath)

