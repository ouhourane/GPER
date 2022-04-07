gper <- function(x, y, group = c(1:nvars),tau=0.5,method=c("GLasso","GMcp","GScad","aGMcp","aGScad"), nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.05, 0.001), 
                lambda = NULL, pf = sqrt(bs), dfmax = as.integer(max(group))+1, 
                pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08, 
                gamm = ifelse(((method == "GMcp")), 3, 4) ,intercept=TRUE) {
  #################################################################################
  this.call <- match.call()
  #\tDesign matrix setup, error checking
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
    stop("The response y must be numeric. Factors must be converted to numeric")
  
  #################################################################################
  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else if (length(group) != nvars) 
    stop("group length does not match the number of predictors in x")
  
  bn <- as.integer(max(group))
  bs <- as.integer(as.numeric(table(group)))
  
  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) 
    stop("Groups must be consecutively numbered 1,2,3,...")
  
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
  for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
  ##for (g in 1:bn) gamma[g] <- 2*max(tau,1-tau)/c
  #################################################################################
  #parameter setup
  if (tau < 0) 
    stop("tau must be non-negtive")
  
  if (length(pf) != bn) 
    stop("The size of group-lasso penalty factor must be same as the number of groups")
  maxit <- as.integer(maxit)
  pf <- as.double(pf)
  gamm <- as.double(gamm)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  r <- y
  #################################################################################
  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1")
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
  intr = 0
  if (intercept) intr <- 1
  intr = as.integer(intr)
  #################################################################################
  # call Fortran core
  gamma <- 2*max(tau,1-tau)*gamma / nobs
  gamma <- as.double(gamma)

    if (method == "aGMcp") {
    fit0 <- .Fortran("asqr1mcp",gamm,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                     b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    
  } 
  if (method == "aGScad") {
    fit0 <- .Fortran("asqr1scad",gamm,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                     b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
  }

  if (method == "GMcp") {
    fit0 <- .Fortran("sqr1mcp",gamm,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                     b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    
  } 
  if (method == "GScad") {
    fit0 <- .Fortran("sqr1scad",gamm,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                     b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
  }
  if (method == "GLasso"){
    fit0 <- .Fortran("sqr1lasso",tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                     as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                     b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                     nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
  }
  #################################################################################
  # output
  fit <- getoutputBeta(fit0, maxit, pmax, nvars, vnames)
  fit <- c(fit, list(npasses = fit0$npass, jerr = fit0$jerr, group = group))
  class(fit) <- c("gperpath")
  if (is.null(lambda)) 
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  class(fit) <- c("gper", class(fit))
  fit
}


