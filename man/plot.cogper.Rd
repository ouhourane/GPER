\name{plot.cogper}
\alias{plot.cogper}
\title{Plot coefficients from a cogper object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted cogper object. This function is modified based on the \code{plot} method in the \code{glmnet} package.
}
\usage{
\method{plot}{cogper}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{cogper}} model}
		\item{xvar}{what is on the x-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. Otherwise, plot the curves with gray colors. Default is \code{FALSE}.}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Otherwise, do not put labels. Default is \code{FALSE}.}
		\item{\dots}{other graphical parameters to plot.}
}
\details{
Two coefficient profile plots are produced, one for the mean coefficients and the other for the scale coefficients.
}
\author{
Yuwen Gu and Hui Zou\cr
Maintainer: Yuwen Gu  <guxxx192@umn.edu>
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
m2 <- cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)
#plot(m2)
}

\seealso{\code{\link{plot.cv.cogper}}}

\keyword{models}
\keyword{regression}
