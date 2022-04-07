\name{plot.cv.gper}
\alias{plot.cv.gper}
\title{Plot the cross-validated curve produced by cv.gper}
\description{
Plots the cross-validated curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv.glmnet} function from the \code{glmnet} package.}

\usage{
\method{plot}{cv.gper}(x, sign.lambda = 1, ...)
}

\arguments{
	\item{x}{fitted \code{\link{cv.gper}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}

\details{A plot is produced.}

\author{Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.90
pf <- abs(rnorm(p))
lambda2 <- 1
m1.cv <- cv.gper(y = y, x = x, tau = tau, eps = 1e-8)
plot(m1.cv)
}

\seealso{\code{\link{cv.gper}}}

\keyword{models}
\keyword{regression}
