\name{nlqmmControl}
\alias{nlqmmControl}

\title{
Control parameters for nlqmm estimation
}
\description{
A list of parameters for controlling the fitting process.
}

\usage{
nlqmmControl(method = "Nelder-Mead", LL_tol = 1e-5, beta = 0.5,
	max_iter = 500,	analytic = FALSE, REoptimizer = "nlm",
	REcontrol = list(),	initialize = "nlme", verbose = FALSE)
}

\arguments{
	\item{method}{
		character vector that specifies the optimization algorithm in \code{\link{optim}} to maximize the log-likelihood. Only \code{"Nelder-Mead"} and \code{"BFGS"} allowed.
	}

	\item{LL_tol}{
		tolerance expressed as absolute change of the log-likelihood.
	}

	\item{beta}{
		decreasing step factor for smoothing parameter \code{omega}. See details.
	}

	\item{max_iter}{
		maximum number of iterations.
	}
	\item{analytic}{
		logical flag. If \code{TRUE}, the gradient and Hessian for computing the modal random effects with \code{\link{nlm}} are calculated analytically. Otherwise, they are computed with numerical differentiation (default).
	}
	\item{REoptimizer}{
		optimizer for the modal random effects. The options are \code{"optim"} for BFGS in \code{\link{optim}} or \code{"nlm"} (default) for Newton-type minimization in \code{\link{nlm}}.
	}
	\item{REcontrol}{
		a list of arguments to be passed to the optimizer for the modal random effects. See arguments in \code{\link{optim}} and \code{\link{nlm}}.
	}
	\item{initialize}{
		character specifying what algorithm should be used to initialize all \code{nlqmm} parameters. The options are \code{"nls"} for \code{\link{nls}}, \code{"nlrq"} for \code{\link{nlrq}}, and \code{"nlme"} (default) for \code{\link{nlme}}. If \code{NULL}, naive estimates are used as initial values.
	}
	\item{verbose}{
		logical flag.
	}

}

\details{
The parameters are initialized using one of \code{"nls"}, \code{"nlrq"}, or \code{"nlme"}. These algorithms are started with the values given via the argument \code{start} in \code{\link{nlqmm}}. Only \code{"nlme"} provides estimates of variance-covariance parameters and random effects. Therefore, when using the other algorithms, these parameters are initialized with a naive estimate.

The parameter \code{omega} controls the quadratic approximation of the absolute deviation function at the kink 0 (Chen, 2007). In \code{\link{nlqmm}}, the starting value for \code{omega} is determined automatically and is not (currently) under the user's control. At each iteration, \code{omega} is decreased by a factor \code{beta}. The smaller \code{omega} at convergence, the smaller the approximation error. See details of the algorithm in Geraci (2017).
}

\value{
a list of control parameters.
}

\references{
Chen C. (2007). A finite smoothing algorithm for quantile regression. Journal of Computational and Graphical Statistics, 16(1), 136-164.

Geraci M (2019). Modelling and estimation of nonlinear quantile regression with clustered data. Computational Statistics and Data Analysis, 136, 30-46.
}

\author{
Marco Geraci}

\seealso{
\code{\link{nlqmm}}
}
\keyword{fitting}
\keyword{control}
