\name{plot.cv.cogper}
\alias{plot.cv.cogper}
\title{Plot the cross-validated curve produced by cv.cogper}
\description{
Plots the cross-validated curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv.glmnet} function from the \code{glmnet} package.}

\usage{
\method{plot}{cv.cogper}(x, sign.lambda = 1, ...)
}

\arguments{
	\item{x}{fitted \code{\link{cv.cogper}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}

\details{A plot is produced.}

\author{
Maintainer: Mohamed Ouhourane <mohamed.ouhourane@gmail.com>
}

\examples{
set.seed(1)
n <- 100
p <- 400
group <- c(1:p)
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
tau <- 0.30
w <- 2.0
lambda2 <- 1
m2.cv <- cv.cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)
plot(m2.cv)
}

\seealso{\code{\link{cv.cogper}}}

\keyword{models}
\keyword{regression}
