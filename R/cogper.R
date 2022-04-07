#w=1; lambda=NULL; method="GLasso"; nlambda = 100; lambda.factor = ifelse(nobs < nvars, 0.05, 0.001); pfscale = sqrt(bs);  pfmean = sqrt(bs); pfL1 = rep(1,nvars); dfmax = as.integer(max(group)) + 1; pmax = min(dfmax * 1.2, as.integer(max(group))); eps = 1e-08; maxit = 3e+08;gamm = 3; tau=0.5; intercept=TRUE

cogper <- function (x, y, w=1, lambda=NULL, group = c(1:nvars), method=c("GLasso","GMcp","GScad","aGMcp","aGScad"), nlambda = 100, 
                    lambda.factor = ifelse(nobs < nvars, 0.05, 0.001), 
                   pf.scale = sqrt(bs),  pf.mean = sqrt(bs), dfmax = as.integer(max(group)) + 1, 
                  pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08, 
                  gamm = ifelse(method == "GMcp", 3, 4), tau=0.5,intercept=TRUE) {
  #################################################################################
  
  this.call <- match.call()
  method <- match.arg(method)
  if (!is.matrix(x)) 
    stop("x has to be a matrix")
  
  if (any(is.na(x))) 
    stop("Missing values in x not allowed!")
  y <- drop(y)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames)) 
    vnames <- paste("V", seq(nvars), sep = "")
  
  if (length(y) != nobs) 
    stop("x and y have different number of rows")
  
  if (!is.numeric(y)) 
    stop("The response y wst be numeric. Factors wst be converted to numeric")
  
  #################################################################################
  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else if (length(group) != nvars) 
    stop("group length does not match the number of predictors in x")
  
  bn <- as.integer(max(group))
  bs <- as.integer(as.numeric(table(group)))
  
  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) 
    stop("Groups wst be consecutively numbered 1,2,3,...")
  
  ix <- rep(NA, bn)
  iy <- rep(NA, bn)
  j <- 1
  for (g in 1:bn) {
    ix[g] <- j
    iy[g] <- j + bs[g] - 1
    j <- j + bs[g]
  }
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)
  #################################################################################
  #  get upper bound
  gamma <- rep(NA, bn)
  for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)/nobs
  #################################################################################
  #parameter setup
  if (tau < 0) 
    stop("tau wst be non-negtive")
  tau <- as.double(tau)
  if (length(pf.mean) != bn) 
    stop("The size of group-lasso penalty factor wst be same as the number of groups")
  if (length(pf.scale) != bn) 
    stop("The size of group-lasso penalty factor wst be same as the number of groups")
  maxit <- as.integer(maxit)
  pfmean <- as.double(pf.mean)
  pfscale <- as.double(pf.scale)
  gamm <- as.double(gamm)
  w <- as.double(w)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ########################################
  #lambda setup
  nlam <- as.integer(nlambda)
  
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambdai.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)
  #########################################
  # call Fortran core
  gamma <- as.double(gamma)
  
  if (method == "GMcp") {
  fit0 <- .Fortran("CserMcpUnifier",gamm,w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), theta = double(nvars * nlam), 
                  idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  } 
  if (method == "aGMcp") {
  fit0 <- .Fortran("aCserMcpUnifier",gamm,w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), theta = double(nvars * nlam), 
                  idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  } 
  if (method == "GScad") {
  fit0 <- .Fortran("CserScadUnifier",gamm,w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), theta = double(nvars * nlam), 
                  idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  } 
  if (method == "aGScad") {
  fit0 <- .Fortran("aCserScadUnifier",gamm,w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), theta = double(nvars * nlam), 
                  idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))
  } 
  if (method == "GLasso") {
  fit0 <- .Fortran("CserLassoUnifier",w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                  as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                  b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), theta = double(nvars * nlam), 
                  idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),ntheta = integer(nlam), 
                  alam = double(nlam), npass = integer(1), jerr = integer(1))

  }
  #################################################################################
  fit <- getoutputThetaBeta(fit0, maxit, pmax, nvars, vnames)
  fit <- c(fit, list(npasses = fit0$npass, jerr = fit0$jerr, group = group))
  class(fit) <- c("cogperpath")
  if (is.null(lambda)) 
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  class(fit) <- c("cogper", class(fit))
  fit
} 


