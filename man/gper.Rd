\name{gper}
\alias{gper}
\title{Regularization paths for the sparse asymmetric least squares (SALES) regression (or the sparse expectile regression)}

\description{Fits regularization paths for the Lasso or elastic net penalized asymmetric least squares regression at a sequence of regularization parameters.
}

\usage{
gper(x, y, nlambda = 100L, method = "er",
      lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04),
      lambda = NULL, lambda2 = 0, pf = rep(1, nvars),
      pf2 = rep(1, nvars), exclude, dfmax = nvars + 1,
      pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
      intercept = TRUE, eps = 1e-08, maxit = 1000000L,
      tau = 0.5)
}

\arguments{
		\item{x}{matrix of predictors, of dimension (nobs * nvars); each row is an observation.}

		\item{y}{response variable.}

		\item{nlambda}{the number of \code{lambda} values (default is 100).}

		\item{method}{a character string specifying the loss function to use. only \code{er} is available now.}

		\item{lambda.factor}{The factor for getting the minimal lambda in the \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p} (the number of predictors). If \eqn{N < p}, the default is \code{0.01}. If \eqn{N > p}, the default is \code{0.0001}, closer to zero.  A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is a user-supplied \code{lambda} sequence.}

		\item{lambda}{a user-supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. It is better to supply, if necessary, a decreasing sequence of \code{lambda} values than a single (small) value. The program will ensure that the user-supplied \code{lambda} sequence is sorted in decreasing order before fitting the model.}
		\item{lambda2}{regularization parameter \code{lambda2} for the quadratic penalty of the coefficients.}

		\item{pf}{L1 penalty factor of length \eqn{p} used for the adaptive LASSO or adaptive elastic net. Separate L1 penalty weights can be applied to each coefficient to allow different L1 shrinkage. Can be 0 for some variables, which imposes no shrinkage, and results in that variable always be included in the model. Default is 1 for all variables (and implicitly infinity for variables listed in \code{exclude}).}

		\item{pf2}{L2 penalty factor of length \eqn{p} used for adaptive elastic net. Separate L2 penalty weights can be applied to each coefficient to allow different L2 shrinkage. Can be 0 for some variables, which imposes no shrinkage. Default is 1 for all variables.}

		\item{exclude}{indices of variables to be excluded from the model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{the maximum number of variables allowed in the model. Useful for very large \eqn{p} when a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{the maximum number of coefficients allowed ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters the model through the path, it will be counted only once. Default is \code{min(dfmax*1.2, p)}.}

		\item{standardize}{logical flag for variable standardization, prior to fitting the model sequence. The coefficients are always returned to the original scale. Default is \code{TRUE}.}

		\item{intercept}{Should intercept(s) be fitted (default is \code{TRUE}) or set to zero (\code{FALSE})?}

		\item{eps}{convergence threshold for coordinate descent. Each inner coordinate descent loop continues until the maximum change in any coefficient is less than \code{eps}. Defaults value is \code{1e-8}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda values. Default is 1e7. If the algorithm does not converge, consider increasing \code{maxit}.}

		\item{tau}{the parameter \eqn{\tau} in the ALS regression model. The value must be in (0,1). Default is 0.5.}
}

\details{
Note that the objective function in \code{gper} is \deqn{\Psi_{\tau}(y-X\beta))/N + \lambda_{1}*|\beta| + 0.5*\lambda_{2}*|\beta|^2,}{\Psi\tau(y-X\beta))/N + \lambda1*|\beta| + 0.5*\lambda2*||\beta||^2,} where \eqn{\Psi_{\tau}}{\Psi\tau} the asymmetric squared error loss and the penalty is a combination of weighted L1 and L2 terms.

For faster computation, if the algorithm is not converging or running slow, consider increasing \code{eps}, decreasing \code{nlambda}, or increasing \code{lambda.factor} before increasing \code{maxit}.
}


\value{
An object with S3 class \code{\link{gper}}.
		\item{call}{the call that produced this object}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of \code{lambda}.}
		\item{dim}{dimension of coefficient matrix}
		\item{npasses}{total number of iterations summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{
Maintainer: Mohamed Ouhourane <mohamed.ouhourane@gmail.com>
}

\references{
Gu, Y. and Zou, H. (Preprint), "High-dimensional Generalizations of Asymmetric Least Squares Regression and Their Applications". \emph{Annals of Statistics}.\cr
}

\seealso{\code{\link{plot.gper}}}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.90
m1 <- gper(y = y, x = x, tau = tau, eps = 1e-8)

}

\keyword{models}
\keyword{regression}
