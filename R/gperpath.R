gperpath <- function(x, y, group, method, gamma, gamm,ix, iy, bn, bs, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                    pf, maxit, tau, nobs, nvars, vnames) {
    #################################################################################
    # data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    if (tau <= 0 || tau >= 1) stop("tau must be in (0,1)")
	  tau <- as.double(tau)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("sqr1lasso",tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
                                as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
                                b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
                                nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1),package="GPER")
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("gperpath")
    outlist
} 
