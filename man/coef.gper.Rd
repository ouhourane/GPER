\name{coef.gper}
\alias{coef.gper}
\alias{coef.gperpath}
\title{Get coefficients from an gper object}
\description{
Computes the coefficients or returns a list of the indices of the nonzero coefficients at the requested values for \code{lambda} from a fitted gper object.
}
\usage{
\method{coef}{gper}(object, s = NULL, type = c("coefficients", "nonzero"), ...)
}
\arguments{
		\item{object}{fitted \code{\link{gper}} model object.}

		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are to be made. Default is the entire sequence used to
		create the model.}

		\item{type}{type \code{"coefficients"} computes coefficients at the requested values for \code{s}. Type \code{"nonzero"} returns a list of the indices of nonzero coefficients for each value of \code{s}. Default is \code{"coefficients"}.}

		\item{\dots}{not used. Other arguments to predict. }
}

\details{
\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{coef} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of coefficients from both left and right \code{lambda} indices.
}
\value{The object returned depends on type.}


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
as.vector(coef(m1, s = m1$lambda[5]))
}

\seealso{\code{\link{gper}}, \code{\link{predict.gper}}}

\keyword{models}
\keyword{regression}
