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
Geraci M (2019). Modelling and estimation of nonlinear quantile regression with clustered data. Computational Statistics and Data Analysis, 136, 30-46.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{nlqmm}}
}

\keyword{random effects}
