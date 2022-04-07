\name{predict.gper}
\alias{predict.gper}
\alias{predict.gperpath}
\title{Make predictions from an gper object}

\description{
Similar to other predict methods, this functions predicts fitted values from a fitted gper object.
}

\usage{
\method{predict}{gper}(object, newx, s = NULL, type = "response", ...)
}

\arguments{
		\item{object}{fitted \code{\link{gper}} model object.}
    
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

\seealso{\code{\link{gper}}, \code{\link{coef.gper}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
