\name{predict.cv.gper}
\alias{predict.cv.gper}
\title{make predictions from a cv.gper object}
\description{
This function makes predictions from a cross-validated gper model,
using the fitted \code{cv.gper} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{predict}{cv.gper}(object, newx, s=c("lambda.1se", "lambda.min"), ...)}

\arguments{
		\item{object}{fitted \code{\link{cv.gper}} object.}

		\item{newx}{matrix of new values for \code{x} at which predictions are to be made. Must be a matrix. See documentation for \code{predict.gper}.}

		\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are to be made. Default is the value \code{s = "lambda.1se"} stored on the CV object. Alternatively \code{s = "lambda.min"} can be used. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}

		\item{\dots}{not used. Other arguments to predict.}} \details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{gper}} objects.}

\author{Maintainer: Mohamed Ouhourane <mohamed.ouhourane@gmail.com>
}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.90
pf <- abs(rnorm(p))
m1.cv <- cv.gper(y = y, x = x, tau = tau, eps = 1e-8)
as.vector(predict(m1.cv, newx = x, s = "lambda.min"))
}

\seealso{\code{\link{cv.gper}}, \code{\link{coef.cv.gper}}}

\keyword{models}
\keyword{regression}
