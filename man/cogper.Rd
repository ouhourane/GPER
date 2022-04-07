\name{cogper}
\alias{cogper}
\title{Regularization paths for the coupled sparse asymmetric least squares (COSALES) regression (or the coupled sparse expectile regression)}
\description{Fits regularization paths for coupled sparse asymmetric least squares regression at a sequence of regularization parameters.}
\usage{
cogper(x, y, w = 1.0, nlambda = 100L, method = "cper", 
        lambda.factor = ifelse(2 * nobs < nvars, 1e-02, 1e-04), 
        lambda = NULL, lambda2 = 0, pf.mean = rep(1, nvars), 
        pf2.mean = rep(1, nvars), pf.scale = rep(1, nvars),
        pf2.scale = rep(1, nvars), exclude, dfmax = nvars + 1, 
        pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
        intercept = TRUE, eps = 1e-08, maxit = 1000000L, 
        tau = 0.80)
}
\arguments{
		\item{x}{matrix of predictors, of dimension (nobs * nvars); each row is an observation.}

		\item{y}{response variable.}
    
    \item{w}{weight applied to the asymmetric squared error loss of the mean part. See details. Default is 1.0.}

		\item{nlambda}{the number of \code{lambda} values (default is 100).}

		\item{method}{a character string specifying the loss function to use. Only \code{cper} is available now.}

		\item{lambda.factor}{The factor for getting the minimal lambda in the \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)} with \code{max(lambda)} being the smallest value of \code{lambda} for which all coefficients are zero. The default value depends on the relationship between \eqn{N} (the number of observations) and \eqn{p} (the number of predictors). If \eqn{N < p}, the default is \code{0.01}. If \eqn{N > p}, the default is \code{0.0001}, closer to zero.  A very small value of \code{lambda.factor} will lead to a saturated fit. The argument takes no effect if there is a user-supplied \code{lambda} sequence.} 

		\item{lambda}{a user-supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}. It is better to supply, if necessary, a decreasing sequence of \code{lambda} values than a single (small) value. The program will ensure that the user-supplied \code{lambda} sequence is sorted in decreasing order.}
    
		\item{lambda2}{regularization parameter \code{lambda2} for the quadratic penalty of the coefficients. Default is 0, meaning no L2 penalization.}

		\item{pf.mean, pf.scale}{L1 penalty factor of length \eqn{p} used for adaptive LASSO or adaptive elastic net. Separate L1 penalty weights can be applied to each mean or scale coefficient to allow different L1 shrinkage. Can be 0 for some variables, which imposes no shrinkage and results in that variable being always included in the model. Default is 1 for all variables (and implicitly infinity for variables listed in \code{exclude}).}

		\item{pf2.mean, pf2.scale}{L2 penalty factor of length \eqn{p}{p} used for adaptive elastic net. Separate L2 penalty weights can be applied to each mean or scale coefficient to allow different L2 shrinkage. Can be 0 for some variables, which imposes no shrinkage. Default is 1 for all variables.}

		\item{exclude}{indices of variables to be excluded from the model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{limit the maximum number of variables in the model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters the model through the path, it will be counted only once. Default is \code{min(dfmax*1.2, p)}.}

		\item{standardize}{logical flag for variable standardization, prior to fitting the model sequence. The coefficients are always returned to the original scale. Default is \code{TRUE}.}

		\item{intercept}{Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).}

		\item{eps}{convergence threshold for coordinate descent. Each inner coordinate descent loop continues until the maximum change in any coefficient is less than \code{eps}. Defaults value is \code{1e-8}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda values. Default is 1e7. If the algorithm does not converge, consider increasing \code{maxit}.}

		\item{tau}{the parameter \code{tau} in the coupled ALS regression model. The value must be in (0,1) and cannot be 0.5. Default is 0.8.}
}

\details{
Note that the objective function in \code{cogper} is \deqn{w*\Psi(y-X\beta,0.5)/N + \Psi(y-X\beta-X\theta,\tau)/N + \lambda_1*|\beta| + 0.5*\lambda_2*||\beta||^2 + \mu_1*|\theta| + 0.5*\mu_2*||\theta||^2,}{w*\Psi(y-X\beta,0.5)/N + \Psi(y-X\beta-X\theta,\tau)/N + \lambda1*|\beta| + 0.5*\lambda2*||\beta||^2 + \mu1*|\theta| + 0.5*\mu2*||\theta||^2,} where \eqn{\Psi(u,\tau)=|\tau-I(u<0)|*u^2} denotes the asymmetric squared error loss and the penalty is a combination of L1 and L2 terms for both the mean and scale coefficients.

For faster computation, if the algorithm is not converging or running slow, consider increasing \code{eps}, decreasing \code{nlambda}, or increasing \code{lambda.factor} before increasing \code{maxit}.
}

\value{
An object with S3 class \code{\link{cogper}}.
		\item{call}{the call that produced this object.}
		\item{b0, t0}{intercept sequences both of length \code{length(lambda)} for the mean and scale respectively.}
		\item{beta, theta}{\code{p*length(lambda)} matrices of coefficients for the mean and scale respectively, stored as sparse matrices (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package). To convert them into normal R matrices, use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df.beta, df.theta}{the number of nonzero mean and scale coefficients respectively for each value of \code{lambda}.}
		\item{dim}{dimensions of coefficient matrices.}
		\item{npasses}{total number of iterations summed over all lambda values.}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{
Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu <guxxx192@umn.edu>
}

\references{
Gu, Y. and Zou, H. (Preprint), "High-dimensional Generalizations of Asymmetric Least Squares Regression and Their Applications". \emph{Annals of Statistics}.\cr
}

\seealso{\code{\link{plot.cogper}}}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.30
w <- 2.0
m2 <- cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)
}

\keyword{models}
\keyword{regression}
