\name{predict.nlqmm}
\alias{predict.nlqmm}

\title{
Predictions from an \code{nlqmm} Object
}

\description{
The predictions at level 0 correspond to predictions based only on the fixed effects estimates. The predictions at level 1 are obtained by adding the modes of the random effects.
}

\usage{
\method{predict}{nlqmm}(object, level = 0, newdata = NULL, ...)
}

\arguments{
  \item{object}{
an \code{nlqmm} object.
}
  \item{level}{
not used.
}
\item{newdata}{
	not used.
}

  \item{\dots}{
not used.
}
}

\value{
a vector of predictions.
}

\references{
Geraci M (2019). Modelling and estimation of nonlinear quantile regression with clustered data. Computational Statistics and Data Analysis, 136, 30-46.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{nlqmm}}
}

\keyword{prediction}
