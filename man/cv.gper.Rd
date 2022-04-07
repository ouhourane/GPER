\name{cv.gper}
\alias{cv.gper}
\title{Cross-validation for gper}
\description{Does k-fold cross-validation for gper, produces a plot, and returns a value for \code{lambda}. This function is based on the \code{cv} function from the \code{glmnet} package.}

\usage{
cv.gper(x, y, lambda = NULL, pred.loss = "loss", 
         nfolds = 5, foldid, tau = 0.5, ...)
}

\arguments{
		\item{x}{\code{x} matrix as in \code{\link{gper}}.}
		\item{y}{response variable \code{y} as in \code{\link{gper}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is \code{NULL}, and \code{\link{gper}} chooses its own sequence.}
		\item{nfolds}{number of folds. Default value is 5. Although 
    \code{nfolds} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowed is 3.}
		\item{foldid}{an optional vector of values between 1 and \code{nfolds}, identifying what fold each observation is in. If supplied, \code{nfolds} will be supressed.}
		\item{pred.loss}{loss function used to calculate cross-validation error. The only option now is \code{"loss"}, which is the asymmetric squared error loss (ASEL).
		}
		\item{tau}{the asymmetry coefficient \eqn{\tau} used in the asymmetric squared error loss.}
		\item{\dots}{other arguments that can be passed to gper.}
}

\details{The function runs \code{\link{gper}} \code{nfolds}+1 times; the first to get the \code{lambda} sequence, and the remainder to compute the fit with each of the folds removed. The average error and standard deviation over the folds are computed. 
}

\value{an object of class \code{\link{cv.gper}} is returned, which is a list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length \code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting purposes).}
		\item{gper.fit}{a fitted \code{\link{gper}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
}

\author{
Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>
}

\seealso{\code{\link{gper}}}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.90
lambda2 <- 1
m1.cv <- cv.gper(y = y, x = x, tau = tau, eps = 1e-8)
}

\keyword{models}
\keyword{regression}

