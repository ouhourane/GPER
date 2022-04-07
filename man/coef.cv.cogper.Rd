\name{coef.cv.cogper}
\alias{coef.cv.cogper}
\title{Get coefficients from a cv.cogper object}
\description{
This function gets coefficients from a cross-validated cogper model,
using the fitted \code{cv.cogper} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
\method{coef}{cv.cogper}(object, s = c("lambda.1se", "lambda.min"),...)
}

\arguments{
	\item{object}{fitted \code{\link{cv.cogper}} object.}
  
	\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the value \code{s="lambda.1se"} stored on the CV \code{object}, it is the largest value of \code{lambda} such that error is within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be used, it is the optimal value of \code{lambda} that gives minimum cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be used.}
  
	\item{\dots}{not used. Other arguments to predict. } 
}

\details{This function makes it easier to use the results of
cross-validation to get coefficients or make coefficient predictions.}
\value{The object returned depends the \dots argument which is passed on to the \code{\link{predict}} method for \code{\link{cogper}} objects.}

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
tau <- 0.30
pf <- abs(rnorm(p))
pf2 <- abs(rnorm(p))
w <- 2.0
m2.cv <- cv.cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)

#as.vector(coef(m2.cv, s = "lambda.min")$beta)
#as.vector(coef(m2.cv, s = "lambda.min")$theta)
}

\seealso{\code{\link{cv.cogper}}, \code{\link{predict.cv.cogper}}} 

\keyword{models}
\keyword{regression}
