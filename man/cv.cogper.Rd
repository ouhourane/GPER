\name{cv.cogper}
\alias{cv.cogper}
\title{Cross-validation for cogper}
\description{Does k-fold cross-validation for \code{cogper}, produces a plot, and returns a value for \code{lambda}. This function is based on the \code{cv} function from the \code{glmnet} package.}

\usage{
cv.cogper(x, y, w = 1.0, lambda = NULL, pred.loss = "loss",
           nfolds = 5, foldid, tau = 0.8, ...)
}

\arguments{
		\item{x}{\code{x} matrix as in \code{\link{cogper}}.}
		\item{y}{response variable \code{y} as in \code{\link{cogper}}.}
    \item{w}{weight applied to the asymmetric squared error loss of the mean part. Default is 1.0.}
		\item{lambda}{optional user-supplied lambda sequence; default is \code{NULL}, and \code{\link{cogper}} chooses its own sequence.}
		\item{nfolds}{number of folds. Default value is 5. Although \code{nfolds} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowed is 3.}
		\item{foldid}{an optional vector of values between 1 and \code{nfolds}, identifying what fold each observation is in. If supplied, \code{nfolds} will be supressed.}
		\item{pred.loss}{loss function used to calculate cross-validation error. The only option now is \code{"loss"}, which is the asymmetric squared error loss (ASEL).}
		\item{tau}{the asymmetry coefficient \eqn{\tau} used in the asymmetric squared error loss.}
		\item{\dots}{other arguments that can be passed to cogper.}
}

\details{The function runs \code{\link{cogper}} \code{nfolds}+1 times. The first gets the \code{lambda} sequence, and the remainder fits the model with each of the folds removed. The average error and standard deviation over the folds are computed. 
}

\value{an object of class \code{\link{cv.cogper}} is returned, which is a list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length \code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{a list of two components, each representing the number of non-zero coefficients at each \code{lambda} in the mean and scale part.}
		\item{name}{a text string indicating type of measure (for plotting purposes).}
		\item{cogper.fit}{a fitted \code{\link{cogper}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
}

\author{
Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>
}

\seealso{\code{\link{cogper}}}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.30
w <- 2.0
m2.cv <- cv.cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)
}
\keyword{models}
\keyword{regression}

