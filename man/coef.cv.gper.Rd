\name{coef.cv.gper}
\alias{coef.cv.gper}
\title{Get coefficients from a cv.gper object}
\description{
This function gets coefficients from a cross-validated gper model,
using the fitted \code{cv.gper} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{coef}{cv.gper}(object, s = c("lambda.1se", "lambda.min"),...)
}

\arguments{
	\item{object}{fitted \code{\link{cv.gper}} object.}
  
	\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the value \code{s="lambda.1se"} stored on the CV \code{object}, it is the largest value of \code{lambda} such that error is within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be used, it is the optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}
  
	\item{\dots}{not used. Other arguments to predict. } 
}

\details{This function makes it easier to use the results of
cross-validation to get coefficients or make coefficient predictions.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{gper}} objects.}

\author{ 
Yuwen Gu and Hui Zou\cr
Maintainer:  Yuwen Gu <guxxx192@umn.edu>
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
as.vector(coef(m1.cv, s = "lambda.min"))
}

\seealso{\code{\link{cv.gper}}, \code{\link{predict.cv.gper}}} 

\keyword{models}
\keyword{regression}
