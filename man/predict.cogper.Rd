\name{predict.cogper}
\alias{predict.cogper}
\alias{predict.cogperpath}
\title{Make predictions from a cogper object}

\description{
Similar to other predict methods, this function predicts fitted values from a cogper object.
}

\usage{
\method{predict}{cogper}(object, newx, s = NULL, type = "response", ...)
}

\arguments{
		\item{object}{fitted \code{\link{cogper}} model object.}
    
		\item{newx}{matrix of new values for \code{x} at which predictions are to be made. NOTE: \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
    
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are to be made. Default is the entire sequence used to
		create the model.}
    
		\item{type}{
		type of prediction required. Only \code{response} is available. Gives predicted response for regression problems.}
    
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are to be made. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


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
predict(m2, newx = x, s = m2$lambda[50])
}

\seealso{\code{\link{cogper}}, \code{\link{coef.cogper}}}

\keyword{models}
\keyword{regression}
