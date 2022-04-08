cogperpath <- function(x, y,group, method, gamma,gamm, ix, iy, bn, bs, w, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd,
                      pfmean, pfscale, maxit,tau, nobs, nvars, vnames) {
    #################################################################################
    #data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    if (tau <= 0 || tau >= 1 || tau == 0.5)
        stop("tau must be in (0,1)\\{0.5}")
	tau <- as.double(tau)
    if (w <= 0) stop("Weight must be positive.")
    w <- as.double(w)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("CserLassoUnifier",w,tau, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x),
                     as.double(y), pfmean,pfscale, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1),
                     b0 = double(nlam), beta = double(nvars * nlam),t0 = double(nlam), phi = double(nvars * nlam),
                     idx = integer(pmax),idxf = integer(pmax),nbeta = integer(nlam),nphi = integer(nlam),
                     alam = double(nlam), npass = integer(1), jerr = integer(1),PACKAGE = "GPER")
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("cogperpath")
    outlist
}
