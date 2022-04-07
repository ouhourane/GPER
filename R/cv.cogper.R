cv.cogper <- function(x, y, group = c(1:nvars), w = 1.0, lambda = NULL, pred.loss = "loss",
                       nfolds = 5, foldid, tau = 0.8, pf.scale = sqrt(bs),  pf.mean = sqrt(bs), method=c("GLasso","GMcp","GScad","aGMcp","aGScad"), ...) {
    pred.loss <- match.arg(pred.loss)
    method <- match.arg(method)
    y <- drop(y)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    ########################################
    bs <- as.integer(as.numeric(table(group)))
    #######################################
    # Fit the model once to get lambda
    cogper.object <- cogper(x, y, group = group, w = w, lambda = lambda, tau = tau, pf.scale = pf.scale,  pf.mean = pf.mean, method = method, ...)
    lambda <- cogper.object$lambda
    # Obtain active set size
    nz <- coef(cogper.object, type = "nonzero")
    nz[[1]] <- sapply(nz[[1]], length)
    nz[[2]] <- sapply(nz[[2]], length)
    if (missing(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = nobs))  
    } else nfolds <- max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be at least 3; nfolds=10 recommended")
    outlist <- vector("list", length = nfolds)
    # Fit the nfolds models
    for (i in seq(nfolds)) {
      whichfold <- (foldid == i)
      y_sub <- y[!whichfold]
      outlist[[i]] <- cogper(x = x[!whichfold, , drop = FALSE], y = y_sub, group = group, method = method,
                              w = w, lambda = lambda, tau = tau, pf.scale = pf.scale,  pf.mean = pf.mean, ...)
    }
    # Calculate pred.loss and the model fit
    fun <- paste("cv", class(cogper.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, w, tau))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
                cvlower = cvm - cvsd, nzero = as.list(nz), name = cvname,
                cogper.fit = cogper.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.cogper"
    obj
} 
