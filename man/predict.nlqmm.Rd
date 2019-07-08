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
Geraci M (2017). Nonlinear quantile mixed models. ArXiv e-prints, 1712.09981. Available at \url{https://arxiv.org/abs/1712.09981}.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{nlqmm}}
}

\keyword{prediction}
