\name{plot.gper}
\alias{plot.gper}
\title{Plot coefficients from an gper object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted gper object. This function is modified based on the \code{plot} method in the \code{glmnet} package.
}
\usage{
\method{plot}{gper}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{gper}} model}
		\item{xvar}{what is on the x-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. Otherwise, plot the curves with gray colors. Default is \code{FALSE}.}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Otherwise, do not put labels. Default is \code{FALSE}.}
		\item{\dots}{other graphical parameters to plot.}
}
\details{
A coefficient profile plot is produced.
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
tau <- 0.90
m1 <- gper(y = y, x = x, tau = tau, eps = 1e-8)
plot(m1)
}

\seealso{\code{\link{plot.cv.gper}}}

\keyword{models}
\keyword{regression}
