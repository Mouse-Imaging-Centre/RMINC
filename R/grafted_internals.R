#### Source code taken from Lme4 and LmerTest

## LmerTest
rhoInitJSS <-
  function(model){
    rho <- new.env(parent = emptyenv())
    rho$model <- model
    rho$y <- getY(model)
    rho$X <- getX(model)
    chol(rho$XtX <- crossprod(rho$X))
    rho$REML <- getREML(model)
    if (class(model) == "lmerMod") 
      rho$s <- summary(model)
    rho$fixEffs <- lme4::fixef(model)
    rho$sigma <- lme4::sigma(model)
    rho$vlist <- sapply(model@cnms, length)
    rho$thopt <- lme4::getME(model, "theta")
    rho$param <- as.data.frame(VarCorr(model))[, "sdcor"]
    return(rho)
  }

getX <-
  function(model){
    return(getME(model, "X"))
  }

getY <-
  function(model){
    return(getME(model, "y"))
  }

devfun5 <-
  function(fm, reml = TRUE) {
    stopifnot(is(fm, "merMod"))
    np <- length(fm@pp$theta)
    nf <- length(lme4::fixef(fm))
    if (!lme4::isGLMM(fm)) 
      np <- np + 1L
    n <- nrow(fm@pp$V)
    ff <- updateModel(fm, . ~ ., getREML(fm), attr(model.matrix(fm), 
                                                   "contrasts"), devFunOnly.lmerTest.private = TRUE)
    envff <- environment(ff)
    if (isLMM(fm)) {
      ans <- function(thpars) {
        stopifnot(is.numeric(thpars), length(thpars) == np)
        ff(thpars[-np])
        sigsq <- thpars[np]^2
        dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1))/sigsq + 
          n * log(2 * pi * sigsq)
        if (reml) {
          p <- ncol(envff$pp$RX())
          dev <- dev + 2 * determinant(envff$pp$RX())$modulus - 
            p * log(2 * pi * sigsq)
        }
        return(dev)
      }
    }
    attr(ans, "thopt") <- fm@pp$theta
    class(ans) <- "devfun5"
    ans
  }

updateModel <-
  function (model, mf.final, reml.lmerTest.private, l.lmerTest.private.contrast, 
            devFunOnly.lmerTest.private = FALSE){
    if (!mf.final == as.formula(paste(".~."))) {
      inds <- names(l.lmerTest.private.contrast) %in% attr(terms(as.formula(mf.final)), 
                                                           "term.labels")
      l.lmerTest.private.contrast <- l.lmerTest.private.contrast[inds]
    }
    nfit <- update(object = model, formula. = mf.final, REML = reml.lmerTest.private, 
                   contrasts = l.lmerTest.private.contrast, devFunOnly = devFunOnly.lmerTest.private, 
                   evaluate = FALSE)
    env <- environment(formula(model))
    assign("l.lmerTest.private.contrast", l.lmerTest.private.contrast, 
           envir = env)
    assign("reml.lmerTest.private", reml.lmerTest.private, envir = env)
    assign("devFunOnly.lmerTest.private", devFunOnly.lmerTest.private, 
           envir = env)
    nfit <- eval(nfit, envir = env)
    return(nfit)
  }

getREML <-
  function (model) 
  {
    if (inherits(model, "merMod")) 
      return(getME(model, "is_REML"))
  }

myhess <-
  function(fun, x, fx = NULL, delta = 1e-04, ...){
    nx <- length(x)
    fx <- if (!is.null(fx)) 
      fx
    else fun(x, ...)
    H <- array(NA, dim = c(nx, nx))
    for (j in 1:nx) {
      xadd <- xsub <- x
      xadd[j] <- x[j] + delta
      xsub[j] <- x[j] - delta
      H[j, j] <- (fun(xadd, ...) - 2 * fx + fun(xsub, ...))/delta^2
      for (i in 1:nx) {
        if (i >= j) 
          break
        xaa <- xas <- xsa <- xss <- x
        xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
        xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
        xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
        xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
        H[i, j] <- (fun(xaa, ...) - fun(xas, ...) - fun(xsa, 
                                                        ...) + fun(xss, ...))/(4 * delta^2)
      }
    }
    H[lower.tri(H)] <- t(H)[lower.tri(H)]
    H
  }

calculateTtestJSS <-
  function(rho, Lc, nrow.res, ddf = "Satterthwaite"){
    resultTtest <- matrix(0, nrow = nrow.res, ncol = 4)
    colnames(resultTtest) <- c("df", "t value", "p-value", "sqrt.varcor")
    if (ddf == "Kenward-Roger") {
      if (!requireNamespace("pbkrtest", quietly = TRUE)) 
        stop("pbkrtest package required for Kenward-Roger's approximations")
      Va <- pbkrtest::vcovAdj(rho$model)
    }
    else {
      vss <- vcovJSStheta2(rho$model)
    }
    for (i in 1:nrow.res) {
      if (ddf == "Kenward-Roger") {
        L <- Lc[, i]
        .ddf <- pbkrtest::get_ddf_Lb(rho$model, L)
        b.hat <- rho$fixEffs
        Lb.hat <- sum(L * b.hat)
        Va.Lb.hat <- t(L) %*% Va %*% L
        t.stat <- as.numeric(Lb.hat/sqrt(Va.Lb.hat))
        p.value <- 2 * pt(abs(t.stat), df = .ddf, lower.tail = FALSE)
        resultTtest[i, 1] <- .ddf
        resultTtest[i, 2] <- t.stat
        resultTtest[i, 3] <- p.value
        resultTtest[i, 4] <- as.numeric(sqrt(Va.Lb.hat))
      }
      else {
        g <- mygrad(function(x) vss(t(Lc[, i]), x), c(rho$thopt, 
                                                      rho$sigma))
        denom <- t(g) %*% rho$A %*% g
        varcor <- vss(t(Lc[, i]), c(rho$thopt, rho$sigma))
        resultTtest[i, 1] <- 2 * (varcor)^2/denom
        resultTtest[i, 2] <- (Lc[, i] %*% rho$fixEffs)/sqrt(varcor)
        resultTtest[i, 3] <- 2 * (1 - pt(abs(resultTtest[i, 
                                                         2]), df = resultTtest[i, 1]))
        resultTtest[i, 4] <- sqrt(varcor)
      }
    }
    return(resultTtest)
  }

mygrad <-
  function(fun, x, delta = 1e-04, method = c("central", "forward", 
                                             "backward"), ...){
    method <- match.arg(method)
    nx <- length(x)
    if (method %in% c("central", "forward")) {
      Xadd <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) + 
        diag(delta, nx)
      fadd <- apply(Xadd, 1, fun, ...)
    }
    if (method %in% c("central", "backward")) {
      Xsub <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) - 
        diag(delta, nx)
      fsub <- apply(Xsub, 1, fun, ...)
    }
    res <- switch(method, forward = (fadd - fun(x, ...))/delta, 
                  backward = (fun(x, ...) - fsub)/delta, central = (fadd - 
                                                                      fsub)/(2 * delta))
    res
  }

vcovJSStheta2 <- 
  function(fm){
    stopifnot(is(fm, "merMod"))
    np <- length(fm@pp$theta)
    nf <- length(lme4::fixef(fm))
    if (!lme4::isGLMM(fm)) 
      np <- np + 1L
    ff2 <- updateModel(fm, . ~ ., getREML(fm), attr(model.matrix(fm), 
                                                    "contrasts"), devFunOnly.lmerTest.private = TRUE)
    envff2 <- environment(ff2)
    if (isLMM(fm)) {
      ans <- function(Lc, thpars) {
        stopifnot(is.numeric(thpars), length(thpars) == np)
        sigma2 <- thpars[np]^2
        ff2(thpars[-np])
        vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi())
        return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))
      }
    }
    class(ans) <- "vcovJSStheta2"
    ans
  }

### lme4

RHSForm <-
  function(form, as.form = FALSE){
    rhsf <- form[[length(form)]]
    if (as.form) 
      reformulate(deparse(rhsf))
    else rhsf
  }

getStart <-
  function(start, lower, pred, returnVal = c("theta", "all")) {
    returnVal <- match.arg(returnVal)
    theta <- pred$theta
    fixef <- pred$delb
    if (!is.null(start)) {
      if (is.numeric(start)) {
        theta <- start
      }
      else {
        if (!is.list(start)) 
          stop("start must be a list or a numeric vector")
        if (!all(sapply(start, is.numeric))) 
          stop("all elements of start must be numeric")
        if (length((badComp <- setdiff(names(start), c("theta", 
                                                       "fixef")))) > 0) {
          stop("incorrect components in start list: ", 
               badComp)
        }
        if (!is.null(start$theta)) 
          theta <- start$theta
        noFixef <- is.null(start$fixef)
        noBeta <- is.null(start$beta)
        if (!noFixef) {
          fixef <- start$fixef
          if (!noBeta) {
            message("Starting values for fixed effects coefficients", 
                    "specified through both 'fixef' and 'beta',", 
                    "only 'fixef' used")
          }
        }
        else if (!noBeta) {
          fixef <- start$beta
        }
      }
    }
    if (length(theta) != length(pred$theta)) 
      stop("incorrect number of theta components (!=", length(pred$theta), 
           ")")
    if (length(fixef) != length(pred$delb)) 
      stop("incorrect number of fixef components (!=", length(pred$delb), 
           ")")
    if (returnVal == "theta") 
      theta
    else c(theta, fixef)
  }

checkConv <-
  function(derivs, coefs, ctrl, lbound, debug = FALSE){
    if (is.null(derivs)) 
      return(NULL)
    if (anyNA(derivs$gradient)) 
      return(list(code = -5L, messages = gettextf("Gradient contains NAs")))
    ntheta <- length(lbound)
    res <- list()
    ccl <- ctrl[[cstr <- "check.conv.grad"]]
    checkCtrlLevels(cstr, cc <- ccl[["action"]])
    wstr <- NULL
    if (doCheck(cc)) {
      scgrad <- tryCatch(with(derivs, solve(chol(Hessian), 
                                            gradient)), error = function(e) e)
      if (inherits(scgrad, "error")) {
        wstr <- "unable to evaluate scaled gradient"
        res$code <- -1L
      }
      else {
        mingrad <- pmin(abs(scgrad), abs(derivs$gradient))
        maxmingrad <- max(mingrad)
        if (maxmingrad > ccl$tol) {
          w <- which.max(maxmingrad)
          res$code <- -1L
          wstr <- gettextf("Model failed to converge with max|grad| = %g (tol = %g, component %d)", 
                           maxmingrad, ccl$tol, w)
        }
      }
      if (!is.null(wstr)) {
        res$messages <- wstr
        switch(cc, warning = warning(wstr), stop = stop(wstr), 
               stop(gettextf("unknown check level for '%s'", 
                             cstr), domain = NA))
      }
      if (!is.null(ccl$relTol) && (max.rel.grad <- max(abs(derivs$gradient/coefs))) > 
          ccl$relTol) {
        res$code <- -2L
        wstr <- gettextf("Model failed to converge with max|relative grad| = %g (tol = %g)", 
                         max.rel.grad, ccl$relTol)
        res$messages <- wstr
        switch(cc, warning = warning(wstr), stop = stop(wstr), 
               stop(gettextf("unknown check level for '%s'", 
                             cstr), domain = NA))
      }
    }
    ccl <- ctrl[[cstr <- "check.conv.singular"]]
    checkCtrlLevels(cstr, cc <- ccl[["action"]])
    if (doCheck(cc)) {
      bcoefs <- seq(ntheta)[lbound == 0]
      if (any(coefs[bcoefs] < ccl$tol)) {
        wstr <- "singular fit"
        res$messages <- c(res$messages, wstr)
        switch(cc, warning = warning(wstr), stop = stop(wstr), 
               stop(gettextf("unknown check level for '%s'", 
                             cstr), domain = NA))
      }
    }
    ccl <- ctrl[[cstr <- "check.conv.hess"]]
    checkCtrlLevels(cstr, cc <- ccl[["action"]])
    if (doCheck(cc)) {
      if (length(coefs) > ntheta) {
        H.beta <- derivs$Hessian[-seq(ntheta), -seq(ntheta)]
        resHess <- checkHess(H.beta, ccl$tol, "fixed-effect")
        if (any(resHess$code != 0)) {
          res$code <- resHess$code
          res$messages <- c(res$messages, resHess$messages)
          wstr <- paste(resHess$messages, collapse = ";")
          switch(cc, warning = warning(wstr), stop = stop(wstr), 
                 stop(gettextf("unknown check level for '%s'", 
                               cstr), domain = NA))
        }
      }
      resHess <- checkHess(derivs$Hessian, ccl$tol)
      if (any(resHess$code != 0)) {
        res$code <- resHess$code
        res$messages <- c(res$messages, resHess$messages)
        wstr <- paste(resHess$messages, collapse = ";")
        switch(cc, warning = warning(wstr), stop = stop(wstr), 
               stop(gettextf("unknown check level for '%s'", 
                             cstr), domain = NA))
      }
    }
    if (debug && length(res$messages) > 0) {
      print(res$messages)
    }
    res
  }

checkCtrlLevels <-
  function (cstr, val, smallOK = FALSE){
    bvals <- c("stop", "warning", "ignore")
    if (smallOK) 
      bvals <- outer(bvals, c("", "Small"), paste0)
    if (!is.null(val) && !val %in% bvals) 
      stop("invalid control level ", sQuote(val), " in ", cstr, 
           ": valid options are {", paste(sapply(bvals, sQuote), 
                                          collapse = ","), "}")
    invisible(NULL)
  }

doCheck <-
  function (x){
    is.character(x) && !any(x == "ignore")
  }

checkHess <-
  function (H, tol, hesstype = ""){
    res <- list(code = numeric(0), messages = list())
    evd <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, 
                    error = function(e) e)
    if (inherits(evd, "error")) {
      res$code <- -6L
      res$messages <- gettextf("Problem with Hessian check (infinite or missing values?)")
    }
    else {
      negative <- sum(evd < -tol)
      if (negative) {
        res$code <- -3L
        res$messages <- gettextf(paste("Model failed to converge:", 
                                       "degenerate", hesstype, "Hessian with %d negative eigenvalues"), 
                                 negative)
      }
      else {
        zero <- sum(abs(evd) < tol)
        if (zero || inherits(tryCatch(chol(H), error = function(e) e), 
                             "error")) {
          res$code <- -4L
          res$messages <- paste(hesstype, "Hessian is numerically singular: parameters are not uniquely determined")
        }
        else {
          res$cond.H <- max(evd)/min(evd)
          if (max(evd) * tol > 1) {
            res$code <- c(res$code, 2L)
            res$messages <- c(res$messages, paste("Model is nearly unidentifiable: ", 
                                                  "very large eigenvalue", "\n - Rescale variables?", 
                                                  sep = ""))
          }
          if ((min(evd)/max(evd)) < tol) {
            res$code <- c(res$code, 3L)
            if (!5L %in% res$code) {
              res$messages <- c(res$messages, paste("Model is nearly unidentifiable: ", 
                                                    "large eigenvalue ratio", "\n - Rescale variables?", 
                                                    sep = ""))
            }
          }
        }
      }
    }
    if (length(res$code) == 0) 
      res$code <- 0
    res
  }
