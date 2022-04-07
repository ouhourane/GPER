\name{predict.cv.cogper}
\alias{predict.cv.cogper}
\title{make predictions from a cv.cogper object}
\description{
This function makes predictions from a cross-validated cogper model,
using the fitted \code{cv.cogper} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{predict}{cv.cogper}(object, newx, s=c("lambda.1se", "lambda.min"), ...)}

\arguments{
		\item{object}{fitted \code{\link{cv.cogper}} object.}
    
		\item{newx}{matrix of new values for \code{x} at which predictions are to be made. Must be a matrix. See documentation for \code{predict.cogper}.}
    
		\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are to be made. Default is the value \code{s = "lambda.1se"} stored on the CV object. Alternatively \code{s = "lambda.min"} can be used. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}
    
		\item{\dots}{not used. Other arguments to predict.}} \details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{cogper}} objects.}

\author{Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}

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
as.vector(predict(m2.cv, newx = x, s = "lambda.min"))
}

\seealso{\code{\link{cv.cogper}}, \code{\link{coef.cv.cogper}}}

\keyword{models}
\keyword{regression}
