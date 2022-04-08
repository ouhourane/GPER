\name{print.gper}
\alias{print.gper}
\title{Print an gper object}

\description{
Print a summary of the gper path at each step along the path.
}

\usage{
\method{print}{gper}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{fitted \code{\link{gper}} object.}
  \item{digits}{significant digits in the output.}
  \item{\dots}{additional print arguments.}
}

\details{
The call that produced the \code{\link{gper}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}

\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.
}

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
tau <- 0.90
pf <- abs(rnorm(p))
m1 <- gper(y = y, x = x, tau = tau, eps = 1e-8)
print(m1)
}

\keyword{models}
\keyword{regression}

