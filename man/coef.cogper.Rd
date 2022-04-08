\name{coef.cogper}
\alias{coef.cogper}
\alias{coef.cogperpath}
\title{Get coefficients from a cogper object}
\description{
Computes the coefficients or returns a list of the indices of the nonzero coefficients at the requested values for \code{lambda} from a fitted cogper object.
}
\usage{
\method{coef}{cogper}(object, s = NULL, type = c("coefficients", "nonzero"), ...)
}
\arguments{
		\item{object}{fitted \code{\link{cogper}} model object.}

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
tau <- 0.30
pf <- abs(rnorm(p))
pf2 <- abs(rnorm(p))
w <- 2.0
m2 <- cogper(y = y, x = x, w = w, tau = tau, eps = 1e-8)
}

\seealso{\code{\link{cogper}}, \code{\link{predict.cogper}}}

\keyword{models}
\keyword{regression}
