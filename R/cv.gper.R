cv.gper <- function(x, y, lambda = NULL,group = c(1:nvars), pred.loss = "loss", 
            nfolds = 5, foldid, tau = 0.5, method=c("GLasso","GMcp","GScad","aGMcp","aGScad"), ...) {
    pred.loss <- match.arg(pred.loss)
    method <- match.arg(method)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    gper.object <- gper(x, y, lambda = lambda, group = group, tau = tau, method = method, ...)
    lambda <- gper.object$lambda
    # predict -> coef
    nz <- sapply(coef(gper.object, type = "nonzero"), length)
    if (missing(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = nobs))  
    } else nfolds <- max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be at least 3; nfolds=10 recommended")
    outlist <- vector("list", length = nfolds)
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
      whichfold <- foldid == i
      y_sub <- y[!whichfold]
      outlist[[i]] <- gper(x = x[!whichfold, , drop = FALSE], 
          y = y_sub, group = group, lambda = lambda, tau = tau, method = method, ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(gper.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
        pred.loss, tau))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlower = cvm - cvsd, nzero = nz, name = cvname, gper.fit = gper.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gper"
    obj
} 
