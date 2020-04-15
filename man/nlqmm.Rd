\name{nlqmm}
\alias{nlqmm}

\title{Fitting Nonlinear Quantile Mixed Models}

\description{
\code{nlqmm} is used to fit nonlinear quantile mixed models based on the asymmetric Laplace distribution.}

\usage{
nlqmm(model, fixed, random, group, covariance = "pdDiag", tau = 0.5,
    data = sys.frame(sys.parent()), subset, weights, na.action = na.fail,
    control = list(), nlmeFit = NULL, start = NULL, gradHess = FALSE,
    fit = TRUE) 

}

\arguments{
  \item{model}{
a nonlinear model formula, with the response on the left of a ~ operator and an expression involving parameters and covariates on the right.
}
  \item{fixed}{
a two-sided linear formula of the form \code{f1+...+fn~x1+...+xm}, or a list of two-sided formulas of the form \code{f1~x1+...+xm}, with possibly different models for different parameters. The \code{f1,...,fn} are the names of parameters included on the right hand side of model and the \code{x1+...+xm} expressions define linear models for these parameters (when the left hand side of the formula contains several parameters, they all are assumed to follow the same linear model, described by the right hand side expression). A \code{1} on the right hand side of the formula(s) indicates a single fixed effects for the corresponding parameter(s).
}
  \item{random}{
a two-sided linear formula of the form \code{r1+...+rn~x1+...+xm}, or a list of two-sided formulas of the form \code{r1~x1+...+xm}, with possibly different models for different parameters. The \code{r1,...,rn} are the names of parameters included on the right hand side of model and the \code{x1+...+xm} expressions define linear models for these parameters (when the left hand side of the formula contains several parameters, they all are assumed to follow the same linear model, described by the right hand side expression). A \code{1} on the right hand side of the formula(s) indicates a single fixed effects for the corresponding parameter(s).
}
  \item{group}{
grouping factor.
}
  \item{covariance}{
variance--covariance matrix of the random effects. Default is \code{pdDiag} (see details).
}
  \item{tau}{
the quantile(s) to be estimated.
}
  \item{data}{
an optional data frame containing the variables named in
   \code{fixed}, \code{random} and \code{group}. By default the variables are taken from the environment from which \code{nlqmm} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
an optional vector specifying likelihood weights (see details).
}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{nlqmm} to print an error message and terminate if there are any
   incomplete observations.
}
\item{control}{
  list of control parameters of the fitting process. See \code{\link{nlqmmControl}}.
}
 \item{nlmeFit}{
a \code{\link{nlme}} object with the same model as specified in the \code{nlqmm} call. If provided, its predicted random effects and estimated variance-covariance matrix will be used as initial values for \code{nlqmm}.
}
 \item{start}{
an optional list of initial (crude) estimates for the fixed effects (see details and \code{\link{nlqmmControl}}). It must be a list. The current implementation does not handle model functions of class selfStart. 
}

\item{gradHess}{
	logical flag. If \code{TRUE} gradient and Hessian of the log-likelihood be calculated numerically with \code{numDeriv}. See details.
}
 \item{fit}{
logical flag. If \code{FALSE} the function returns a list of arguments for fitting.
}
}


\details{
The function computes an estimate on the tau-th quantile function of the response, conditional on the covariates, as specified by the \code{formula} argument, and on random effects, as specified by the \code{random} argument. The quantile predictor is assumed to be nonlinear. The function maximizes the (log)likelihood of the Laplace regression proposed by Geraci (2019) using a Laplacian approximation of the likelihood. The optimization algorithm minimizes the negative smoothed Laplace log--likelihood using \code{\link{optim}}. The available options are the Nelder-Mead and BFGS algorithms which are specified via the argument \code{method} in \code{\link{nlqmmControl}}.

Optional starting values for the fixed effects are provided via \code{start}. Parameters are named internally as \code{b1}, \code{b2}, ... See examples further below. These values are then passed on to the initialization algorithm (see \code{\link{nlqmmControl}}) which provides the initial values to the \code{nlqmm} algorithm.

Different standard types of positive--definite matrices for the random effects can be specified: \code{pdIdent} multiple of an identity; \code{pdCompSymm} compound symmetry structure (constant diagonal and constant off--diagonal elements); \code{pdDiag} diagonal; \code{pdSymm} general positive--definite matrix, with no additional structure.

Weights are assumed to be the same for all observations within the same cluster. Weights must be strictly positive. These represent likelihood weights: the larger the weight for a cluster, the more importance is given to the cluster during estimation.

If \code{gradHess = TRUE}, the gradient and Hessian of the log-likelihood are calculated at the optimum. The results are stored as attributes of \code{logLik} and can be extracted via \code{\link{logLik.nlqmm}}. If calculation fails, a warning will be printed out.
}

\value{
\code{nlqmm} returns an object of \code{\link{class}} \code{nlqmm}.
}

\references{
Chen C. (2007). A finite smoothing algorithm for quantile regression. Journal of Computational and Graphical Statistics, 16(1), 136-164.

Geraci M. (2019). Modelling and estimation of nonlinear quantile regression with clustered data. Computational Statistics and Data Analysis, 136, 30-46.
}

\author{
Marco Geraci
}

\examples{
\dontrun{

# Pharmacokinetics of Indomethacin (package datasets)

Indometh$Subject <- as.factor(as.numeric(Indometh$Subject))
tau <- c(.1,.5,.9)
nt <- length(tau)

# Provide starting values for fixed effects
ss <- list(b1 = 3, b2 = 0.5, b3 = 0.5, b4 = -1)

# Fit bi-exponential mixed models for 3 quantiles 
fit <- list()
ctrl <- nlqmmControl(method = "Nelder-Mead", initialize = "nlrq", verbose = FALSE)
for(i in 1:nt){
	fit[[i]] <- nlqmm(conc ~ A1*exp(-exp(lrc1)*time)+A2*exp(-exp(lrc2)*time),
	data = Indometh, fixed = A1 + lrc1 + A2 + lrc2 ~ 1,
	random = A1 + lrc1 + A2 ~ 1, group = Subject, tau = tau[i],
	covariance = "pdDiag", start = ss, control = ctrl, fit = TRUE)
}
}

\dontrun{
# Growth of soybean plants (package nlme)

data(Soybean, package = "nlme")
tau <- c(0.05, .1, .5, .9, .95)
nt <- length(tau)
ss <- c(20, -8, -3, 1.6, 5.5, 0.1, 54.8, -2, -5, -1, 8, -1, -0.6) # Pinheiro and Bates (2000, p.294)

# Fit logistic mixed models for 5 quantiles 
ctrl <- nlqmmControl(method = "BFGS", initialize = "nlrq", verbose = FALSE,
	beta = 0.2, LL_tol = 1e-05)
names(ss) <- paste0("b", 1:13)
fit <- list()
for(i in 1:nt){
	fit[[i]] <- nlqmm(weight ~ Asym/(1+exp((xmid-Time)/scal)), data = Soybean,
	fixed = list(Asym ~ Year*Variety, xmid ~ Year + Variety, scal ~ Year),
	random = Asym ~ 1, group = Plot, tau = tau[i], covariance = "pdSymm",
	start = as.list(ss), control = ctrl, fit = TRUE)
}

}

}

\note{
Updates/FAQ/news are published here \url{http://marcogeraci.wordpress.com/}. New versions are usually published here \url{https://r-forge.r-project.org/R/?group_id=1396} before going on CRAN.
}

\seealso{
\code{\link{lqmm}}
}

\keyword{quantile regression}
