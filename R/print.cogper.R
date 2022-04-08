print.cogper <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df1 = x$df.beta, Df2 = x$df.phi,
                Lambda = signif(x$lambda, digits)))
}
