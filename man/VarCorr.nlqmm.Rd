\name{VarCorr.nlqmm}
\alias{VarCorr}
\alias{VarCorr.nlqmm}

\title{
Extract Variance-Covariance Matrix
}
\description{
This function extracts the variance-covariance matrix of the random effects from a fitted \code{nlqmm} object.
}
\usage{
\method{VarCorr}{nlqmm}(x, sigma = NULL, ...)
}

\arguments{
  \item{x}{
an object of \code{\link{class}} "nlqmm".
}

 \item{sigma}{
not used.
}

 \item{...}{
not used.
}

}

\details{
\code{nlqmm}: This function returns the variance or the variance-covariance matrix of the random effects. The generic function \code{VarCorr} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}

\references{
Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{nlqmm}}
}

\keyword{covariance}
\keyword{coefficients}
