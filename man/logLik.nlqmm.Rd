\name{logLik.nlqmm}
\alias{logLik.nlqmm}

\title{
Extract Log-Likelihood
}
\description{
\code{logLik.nlqmm} extracts the log-likelihood of a fitted nlqmm.
}
\usage{
\method{logLik}{nlqmm}(object, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "nlqmm".

}
  \item{\dots}{
not used.
}
}

\value{

Returns the loglikelihood of the fitted model. This is a number with at least one attribute, "df" (degrees of freedom), giving the number of (estimated) parameters in the model.

The attributes "grad" and "hessian" will also be given if argument \code{gradHess = TRUE} in \code{nlqmm} main call, and calculation has been successful.

}


\author{
Marco Geraci
}


\seealso{
\code{\link{nlqmm}}
}

