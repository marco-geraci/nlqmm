\name{ranef.nlqmm}
\alias{ranef}
\alias{ranef.nlqmm}

\title{
Predictions from an \code{nlqmm} Object
}

\description{
Random effects.
}

\usage{
\method{ranef}{nlqmm}(object, ...)
}

\arguments{
  \item{object}{
an \code{nlqmm} object.
}
\item{...}{
	not used.
}

}

\details{
The modes of the random effects for the Laplacian approximation of the (log)likelihood as obtained at the last iteration.
}

\value{
a vector of predicted random effects.
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

\keyword{random effects}
