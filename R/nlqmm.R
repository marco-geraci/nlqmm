################################################################################
### Fit a nonlinear quantile mixed models (nlqmm)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

##########################################################################################

#' @export
nlqmm <- function(model, fixed, random, group, covariance = "pdDiag", tau = 0.5, data = sys.frame(sys.parent()), subset, weights, na.action = na.fail, control = list(), nlmeFit = NULL, start = NULL, gradHess = FALSE, fit = TRUE) {

Call <- match.call()
  
if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(length(tau) > 1) stop("One 'tau' at the time")

# check arguments
if(!is.data.frame(data)) stop("`data' must be a data frame")
if (!inherits(model, "formula")){
    stop("'model' must be a formula")
}
if (length(model)!=3){
	stop("model formula must be of the form \"y ~ x\"")
}
if (missing(random)) {
	random <- fixed
}

groupFormula <- asOneSidedFormula(Call[["group"]])
group <- groupFormula[[2]]
namGrp <- as.character(group)

if (!is.list(fixed)) fixed <- list(fixed)
fixed <- do.call(c, lapply(fixed, function(fix.i) {
  if (is.name(fix.i[[2]])) list(fix.i)
  else
	## multiple parameters on left hand side
	eval(parse(text = paste0("list(", paste(paste(all.vars(fix.i[[2]]), deparse (fix.i[[3]]), sep = "~"), collapse = ","), ")")))
}))
fnames <- lapply(fixed, function(fix.i) {
	this <- eval(fix.i)
	if (!inherits(this, "formula"))
	  stop ("'fixed' must be a formula or list of formulae")
	if (length(this) != 3)
	  stop ("formulae in 'fixed' must be of the form \"parameter ~ expr\"")
	if (!is.name(this[[2]]))
	  stop ("formulae in 'fixed' must be of the form \"parameter ~ expr\"")
	as.character(this[[2]])
})
names(fixed) <- fnames

if (!is.list(random)) random <- list(random)
random <- do.call(c, lapply(random, function(ran.i) {
  if (is.name(ran.i[[2]])) list(ran.i)
  else
	## multiple parameters on left hand side
	eval(parse(text = paste0("list(", paste(paste(all.vars(ran.i[[2]]), deparse (ran.i[[3]]), sep = "~"), collapse = ","), ")")))
}))
rnames <- lapply(random, function(ran.i) {
	this <- eval(ran.i)
	if (!inherits(this, "formula"))
	  stop ("'random' must be a formula or list of formulae")
	if (length(this) != 3)
	  stop ("formulae in 'random' must be of the form \"parameter ~ expr\"")
	if (!is.name(this[[2]]))
	  stop ("formulae in 'random' must be of the form \"parameter ~ expr\"")
	as.character(this[[2]])
})
names(random) <- rnames

## all parameter names
pnames <- unique(c(fnames, rnames))

## extract data frame with all necessary information
mfArgs <- list(formula = asOneFormula(model, random, fixed, group, omit = c(pnames, "pi")), data = data, na.action = na.action)
if(!missing(subset)) {
	mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
}
if(!missing(weights)) {
	mfArgs[["weights"]] <- weights
}
mfArgs$drop.unused.levels <- TRUE
dataMix <- do.call("model.frame", mfArgs)
origOrder <- row.names(dataMix)	# preserve the original order

## sort the model.frame by groups
grp <- model.frame(groupFormula, dataMix)

## ordering data by groups
ord <- order(unlist(grp, use.names = FALSE))
grp <- grp[ord,,drop = TRUE]
dataMix <- dataMix[ord, ,drop = FALSE]
revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
ngroups <- length(unique(grp))

## obtaining model matrices and vectors
y <- as.matrix(model.frame(model[1:2], dataMix))

## Likelihood weights
if(!missing(weights)) weights <- model.weights(dataMix)[!duplicated(grp)]
if(!missing(weights) && is.null(weights)) weights <- rep(1,ngroups)
if(missing(weights))  weights <- rep(1,ngroups)

weights <- rep((2/table(grp))*log(weights), table(grp))

mmf <- lapply(fixed, function(x, data) model.matrix(x[c(1,3)], data), data = dataMix)
mmr <- lapply(random, function(x, data) model.matrix(x[c(1,3)], data), data = dataMix)

## define dimensions
cov_name <- covariance
dim_theta <- integer(2)
dim_theta[1] <- sum(sapply(mmf, ncol))
dim_theta[2] <- sum(sapply(mmr, ncol))
dim_theta_z <- theta.z.dim(type = cov_name, n = dim_theta[2])

## Control
if(is.null(names(control))) control <- nlqmmControl()
    else {
    control_default <- nlqmmControl();
    control_names <- intersect(names(control), names(control_default));
    control_default[control_names] <- control[control_names];
    control <- control_default
    }
method <- control$method
if(!method %in% c("Nelder-Mead", "BFGS")) {control$method <- method <- "Nelder-Mead"; cat("Switching to Nelder-Mead optimization \n")}

## Parse model
frmod <- parsef(mmf, mmr, model)
fmod <- attr(frmod, "fmodel")

## New data frame for initial estimates
ndata <- cbind(dataMix,mmf)
nn <- names(ndata)
for(i in 1:length(fnames)) nn <- gsub(paste0(fnames[[i]],"."), "", nn)
names(ndata) <- nn
ndata <- ndata[,!duplicated(nn)]

## Check starting values for fixed effects
if(!is.null(start)){
	if(!is.list(start)) stop("'start' must be a list")
	if(!all(attr(frmod, "coefs")$fixed %in% names(start))){
		cat("Model: ~",  as.character(fmod)[[3]], "\n")
		cat("Parameters:",  attr(frmod, "coefs")$fixed, "\n")
		stop("'start' must include starting values for all fixed effects")
	}
}

## If not provided, obtain crude starting values for fixed effects
if(is.null(start)){
	Dfmod <- deriv3(fmod, attr(frmod, "coefs")$fixed)
	get_initial <- function(b, data, nm, model){
		b <- as.list(b)
		names(b) <- nm
		res <- y - as.numeric(eval(model, envir = c(data, b)))
		sum(abs(res))
	}
	start <- optim(fn = get_initial, par = rep(0, sum(dim_theta[1])), data = ndata, nm = attr(frmod, "coefs")$fixed, model = Dfmod)$par
	start <- as.list(start)
	names(start) <- attr(frmod, "coefs")$fixed
}


## Initialize

init <- control$initialize

if(is.null(init)){
	theta_x <- as.numeric(unlist(start))
	theta_z <- rep(0, dim_theta_z)
	theta_0 <- c(theta_x, theta_z)
	sigma_0 <- 1
	RE_0 <- matrix(0, ngroups, dim_theta[2])
	omega_0 <- sd(y)/2
	init <- "none"
}

if(init == "nls"){
	fit_s <- nls(fmod, data = ndata, start = start)
	theta_x <- coef(fit_s)
	theta_z <- rep(0, dim_theta_z)
	theta_0 <- c(theta_x, theta_z)
	sigma_0 <- sqrt(deviance(fit_s)/length(y))
	RE_0 <- matrix(0, ngroups, dim_theta[2])
	omega_0 <- max(kronecker(abs(fit_s$m$resid()), c(tau, 1- tau), "/"))/2 #/sd(y)
}

if(init == "nlrq"){
	fit_s <- quantreg::nlrq(fmod, data = ndata, start = start, tau = tau)
	theta_x <- coef(fit_s)
	theta_z <- rep(0, dim_theta_z)
	theta_0 <- c(theta_x, theta_z)
	sigma_0 <- sqrt(deviance(fit_s)/length(y))
	RE_0 <- matrix(0, ngroups, dim_theta[2])
	omega_0 <- max(kronecker(abs(fit_s$m$resid()), c(tau, 1- tau), "/"))/2 #/sd(y)
}

if(init == "nlme"){
	if(is.null(nlmeFit)){
		random.nlme <- do.call(cov_name, args = list(form = random))
		start.nlme <- as.numeric(start)
		names(start.nlme) <- names(attr(frmod, "labels")$fixed)
		nlmeFit <- try(nlme(model = model, data = data, fixed = fixed, random = random.nlme, start = start.nlme, groups = groupFormula), silent = TRUE)
		if(!inherits(nlmeFit, "try-error")){
			theta_x <- fixef(nlmeFit)
			theta_z <- coef(nlmeFit[[1]])
			sigma_0 <- sigma(nlmeFit)
			RE_0 <- ranef(nlmeFit)
			omega_0 <- max(kronecker(abs(nlmeFit$residuals[,'fixed']), c(tau, 1- tau), "/"))/2 #/sd(y)

		} else {
			theta_x <- as.numeric(unlist(start))
			theta_z <- rep(0, dim_theta_z)
			sigma_0 <- 1
			RE_0 <- matrix(0, ngroups, dim_theta[2])
			omega_0 <- sd(y)/2
		}
	} else {
		theta_x <- fixef(nlmeFit)
		theta_z <- coef(nlmeFit[[1]])
		sigma_0 <- sigma(nlmeFit)
		RE_0 <- ranef(nlmeFit)
		omega_0 <- max(kronecker(abs(nlmeFit$residuals[,'fixed']), c(tau, 1- tau), "/"))/2 #/sd(y)
	}
	if(any(theta_z < -3.45)){
		warning("Some parameters of the scaled variance-covariance matrix from 'nlme' are very small (< 0.001). Possibly singular matrix.")
		theta_z <- rep(0, dim_theta_z)
	}
	theta_0 <- c(theta_x, theta_z)
}

FIT_ARGS <- InitialPar <- list(theta = theta_0, sigma = sigma_0, ranef = RE_0, omega = omega_0)

## Complete FIT_ARGS list with all necessary arguments

FIT_ARGS <- c(FIT_ARGS, list(model = frmod, data = ndata, y = y, weights = as.numeric(weights), cov_name = cov_name, group = grp, tau = tau, analytic = control$analytic, REoptimizer = control$REoptimizer, REcontrol = control$REcontrol))

if(!fit) return(FIT_ARGS)

## Estimation

iter <- 0
FIT_ARGS$ranef <- do.call(modalRe_nlqmm, args = FIT_ARGS[match(names(formals(modalRe_nlqmm)), names(FIT_ARGS))])
ll_0 <- do.call(loglik_nlqmm, args = FIT_ARGS[match(names(formals(loglik_nlqmm)), names(FIT_ARGS))])
delta <- NA

while(iter < control$max_iter){

if(control$verbose){
	cat("iter:", iter, "\n")
	cat("loglikelihood:", ll_0, "\n")
	cat("delta:", delta, "\n")
	cat("theta:", FIT_ARGS$theta, "\n")
	cat("omega:", FIT_ARGS$omega, "\n")
}
# control = list(reltol = control$LL_tol)
tmp <- do.call(optim, args = c(list(fn = loglik_nlqmm, par = FIT_ARGS$theta, method = method), FIT_ARGS[-c(match(c("theta","analytic","REoptimizer","REcontrol"), names(FIT_ARGS)))]))
FIT_ARGS$theta <- tmp$par
ll_1 <- tmp$value

delta <- abs((ll_1 - ll_0)/(ll_0))

	if(delta < control$LL_tol){
		cat("Convergence reached after ", iter + 1, "iteration(s)", "\n")
		break
	} else {
		ll_0 <- ll_1
		FIT_ARGS$omega <- FIT_ARGS$omega*control$beta
		iter <- iter + 1
	}
}

# update sigma, random effects and loglikelihood
logLik <- do.call(loglik_nlqmm, args = FIT_ARGS[match(names(formals(loglik_nlqmm)), names(FIT_ARGS))])
FIT_ARGS$sigma <- attributes(logLik)$sigma
FIT_ARGS$ranef <- do.call(modalRe_nlqmm, args = FIT_ARGS[match(names(formals(modalRe_nlqmm)), names(FIT_ARGS))])
logLik <- -do.call(loglik_nlqmm, args = FIT_ARGS[match(names(formals(loglik_nlqmm)), names(FIT_ARGS))])

# Gradient and Hessian of loglike is time consuming
if(gradHess){
	attr(logLik, "grad") <- -do.call(numDeriv::grad, args = c(list(func = loglik_nlqmm, x = FIT_ARGS$theta), FIT_ARGS[-c(1)]))
	attr(logLik, "hessian") <- -do.call(numDeriv::hessian, args = c(list(func = loglik_nlqmm, x = FIT_ARGS$theta), FIT_ARGS[-c(1)]))
}

fit <- list()
fit$theta_x <- FIT_ARGS$theta[1:dim_theta[1]];
fit$theta_z <- FIT_ARGS$theta[(dim_theta[1] + 1) : (dim_theta[1] + dim_theta_z)]
fit$sigma <- attributes(logLik)$sigma
fit$ranef <- FIT_ARGS$ranef
fit$logLik <- as.numeric(logLik)
fit$residuals <- attr(logLik, "resid")
fit$fitted <- attr(logLik, "fitted")
fit$omega <- FIT_ARGS$omega
fit$call <- Call
fit$model <- frmod
fit$nn <- attr(frmod, "labels")$fixed
fit$mm <- attr(frmod, "labels")$random
fit$nobs <- length(y)
fit$dim_theta <- dim_theta
fit$dim_theta_z <- dim_theta_z
fit$edf <- fit$dim_theta[1] + fit$dim_theta_z
fit$rdf <- fit$nobs - fit$edf
fit$df <- dim_theta[1] +  dim_theta_z + 1
fit$tau <- tau
fit$y <- y
fit$revOrder <- revOrder
fit$weights <- weights
fit$group <- grp
fit$ngroups <- ngroups
fit$InitialPar <- InitialPar
fit$control <- control
fit$opt <- list(iter = iter, delta = as.numeric(delta))
fit$cov_name <- cov_name
fit$mfArgs <- mfArgs

class(fit) <- "nlqmm"
fit
}

#' @export
nlqmmControl <- function(method = "Nelder-Mead", LL_tol = 1e-5, beta = 0.5,	max_iter = 500,	analytic = FALSE, REoptimizer = "nlm", REcontrol = list(), initialize = "nlme", verbose = FALSE){

if(!method %in% c("Nelder-Mead", "BFGS")) {method <- "Nelder-Mead"; cat("Switching to Nelder-Mead optimization \n")}

if(beta > 1 || beta < 0) stop("Beta must be a decreasing factor in (0,1)")
if(max_iter < 0) stop("Number of iterations cannot be negative")

list(method = method, LL_tol = LL_tol, beta = beta, max_iter = as.integer(max_iter), analytic = analytic, REoptimizer = REoptimizer, REcontrol = REcontrol, initialize = initialize, verbose = verbose)

}

parsef <- function(mmf, mmr, model){

pnames <- unique(c(names(mmf), names(mmr)))
b <- bl <- lapply(mmf, function(x) colnames(x))
u <- ul <- lapply(mmr, function(x) colnames(x))
P <- sapply(mmf, function(x) ncol(x))
Q <- sapply(mmr, function(x) ncol(x))
minp <- c(1, cumsum(P[-length(P)]) + 1)
maxp <- cumsum(P)
minq <- c(1, cumsum(Q[-length(Q)]) + 1)
maxq <- cumsum(Q)
	for(i in 1:length(b)){
		tmp <- paste(paste0("(",b[[i]],")"), paste0("b",seq(minp[i], maxp[i])), sep = "*")
		bl[[i]] <- gsub("\\(\\(Intercept\\)\\)\\*", "", tmp)
	}
	for(i in 1:length(u)){
		tmp <- paste(paste0("(",u[[i]],")"), paste0("u",seq(minq[i], maxq[i])), sep = "*")
		ul[[i]] <- gsub("\\(\\(Intercept\\)\\)\\*", "", tmp)
	}

# rewrite model with fixed and random coefficients
frmod <- as.character(model)[[3]]
	for(i in 1:length(pnames)){
		sel <- pnames[[i]]
		tmp <- paste(c(bl[[match(sel, names(bl))]], ul[[match(sel, names(ul))]]), collapse = " + ")
		frmod <- gsub(sel, paste0("(", tmp, ")"), frmod)
		frmod <- gsub(":", "*", frmod)
	}

# rewrite model with fixed coefficients only, ie u = 0 (to get initial values)
fmod <- as.character(model)[[3]]
	for(i in 1:length(pnames)){
		sel <- pnames[[i]]
		tmp <- paste(bl[[match(sel, names(bl))]], collapse = " + ")
		fmod <- gsub(sel, paste0("(", tmp, ")"), fmod)
		fmod <- gsub(":", "*", fmod)
	}
	
frmod <- formula(paste(as.character(model)[[2]], "~", frmod))
fmod <- formula(paste(as.character(model)[[2]], "~", fmod))

attr(frmod, "coefs") <- list(fixed = paste0("b",1:sum(P)), random = paste0("u",1:sum(Q)))
attr(frmod, "labels") <- list(fixed = b, random = u)
attr(frmod, "fmodel") <- fmod
attr(frmod, "model") <- model
attr(frmod, "mmf") <- mmf
attr(frmod, "mmr") <- mmr

return(frmod)
}

#' @export
hfun_nlqmm <- function(ranef, theta_x, H, model, data, y, tau, omega, analytic){

all_coefs <- as.list(c(theta_x, ranef))
names(all_coefs) <- unlist(attr(model, "coefs"))
Dfrmod <- deriv3(model, attr(model, "coefs")$random)
eta <- eval(Dfrmod, envir = c(data, all_coefs))
res <- y - as.numeric(eta)

ans <- C_hfun_nlqmm(res, ranef, H, length(res), length(ranef), tau, omega)

## Gradient
#if(analytic){
#J <- attr(val, "gradient")
#grad <- -t(J) %*% matrix(2/omega * w * res + bvec) + 2*H %*% ranef
#attr(ans, "gradient") <- grad
#}
return(ans)

}

#' @export
modalRe_nlqmm <- function(ranef, theta, sigma, model, data, y, cov_name, group, tau, omega, analytic, REoptimizer, REcontrol){

eps <- .Machine$double.eps^(2/3)

id <- unique(group)
M <- length(id)
feTrms <- attr(model, "coefs")$fixed
reTrms <- attr(model, "coefs")$random
P <- length(feTrms)
Q <- length(reTrms)

theta_x <- theta[1:P]
theta_z <- theta[-c(1:P)]

Sigma <- switch(cov_name,
		pdIdent = pdIdent(theta_z, nam = reTrms),
		pdDiag = pdDiag(theta_z, nam = reTrms),
		pdSymm = pdSymm(theta_z, nam = reTrms),
		pdCompSymm = pdCompSymm(theta_z, nam = reTrms)
)
Psi <- as.matrix(Sigma)/sigma
Psiinv <- try(solve(Psi), silent = TRUE)
if(inherits(Psiinv, "try-error")) Psiinv <- diag(100, Q, Q)

## nlm args

if(REoptimizer == "nlm"){
	control_default <- list(hessian = TRUE, fscale = 1, print.level = 0, ndigit = 12, gradtol = 1e-3, steptol = 1e-4, iterlim = 100, check.analyticals = FALSE)

	if(is.null(names(REcontrol))){
		REcontrol <- control_default} else {
		control_names <- intersect(names(REcontrol), names(control_default));
		control_default[control_names] <- REcontrol[control_names];
		REcontrol <- control_default
	}

	# Force calculation of Hessian. Remove when analytic gradient is implemented
	REcontrol$hessian <- TRUE
	
	nlmArgs <- c(list(f = hfun_nlqmm, theta_x = theta_x, H = Psiinv, model = model, tau = tau, omega = omega, analytic = analytic), REcontrol)
}

## optim args

if(REoptimizer == "optim"){
	control_default <- list(maxit = 100, reltol = 1e-5, usenumDeriv = TRUE, kkt = FALSE)

	if(is.null(names(REcontrol))){
		REcontrol <- control_default} else {
		control_names <- intersect(names(REcontrol), names(control_default));
		control_default[control_names] <- REcontrol[control_names];
		REcontrol <- control_default
	}

	# Force calculation of Hessian. Remove when analytic gradient is implemented
	REcontrol$hessian <- TRUE

	optimArgs <- c(list(fn = hfun_nlqmm, theta_x = theta_x, H = Psiinv, model = model, tau = tau, omega = omega, analytic = analytic, method = "BFGS", gr = NULL, hessian = TRUE), control = REcontrol)
}

# Loop over clusters
uhat <- hfirst <- matrix(NA, M, Q)
hsec <- rep(NA, M)
for(i in 1:M){
	sel <- group == id[i]
	datas <- subset(data, sel)
	yi <- y[sel]
	ui <- as.numeric(ranef[i,1:Q])
	ui[abs(ui) < eps] <- eps

	if(REoptimizer == "nlm"){
		nlmArgs$p <- ui
		nlmArgs$data <- datas
		nlmArgs$y <- yi
		nlmArgs$typsize <- rep(1, length(ui))
		nlmArgs$stepmax <- max(1000 * sqrt(sum((ui/rep(1, length(ui)))^2)), 1000)
		fit <- do.call(nlm, nlmArgs)

		msg <- switch(as.character(fit$code),
			"1" = "relative gradient is close to zero, current iterate is probably solution",
			"2" = "successive iterates within tolerance, current iterate is probably solution",
			"3" = "last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small",
			"4" = "iteration limit exceeded",
			"5" = "maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small"
		)
	}

	if(REoptimizer == "optim"){
		optimArgs$par <- ui
		optimArgs$data <- datas
		optimArgs$y <- yi
		fit <- do.call(optimx::optimx, optimArgs)
		fit$estimate <- as.numeric(unlist(fit[1:length(ui)]))
		sel <- match("ngatend", dimnames(attr(fit, "details"))[[2]])
		fit$gradient <- attr(fit, "details")[[sel]]
		sel <- match("nhatend", dimnames(attr(fit, "details"))[[2]])
		fit$hessian <- attr(fit, "details")[[sel]]

		msg <- switch(as.character(fit$convergence),
			"0" = "successful completion",
			"1" = "iteration limit maxit had been reached",
			"10" = "degeneracy of the Nelder-Mead simplex",
			"51" = "warning from the 'L-BFGS-B' method",
			"52" = "error from the 'L-BFGS-B' method"
		)

	}

	uhat[i,] <- fit$estimate
	hfirst[i,] <- fit$gradient
	hsec[i] <- as.numeric(determinant(fit$hessian, logarithm = TRUE)$modulus)
}

attr(uhat, "hfirst") <- hfirst
attr(uhat, "hsec") <- hsec
return(uhat)

}

#' @export
loglik_nlqmm <- function(theta, sigma, ranef, model, data, y, weights = weights, cov_name, group, tau, omega = 0.001) {
# Argument ranef is an M x Q matrix of random effects modes,
# gradient, and hessian of h function calculated with modalRe_nlqmm

eps <- .Machine$double.eps^(2/3)

id <- unique(group)
M <- length(id)
z <- model.matrix(~group - 1)

feTrms <- attr(model, "coefs")$fixed
reTrms <- attr(model, "coefs")$random
P <- length(feTrms)
Q <- length(reTrms)
N <- length(y)
uhat <- as.matrix(ranef)
hsec <- attr(ranef, "hsec")
logDetH <- sum(hsec)

theta_x <- theta[1:P]
theta_z <- theta[-c(1:P)]

Sigma <- switch(cov_name,
		pdIdent = pdIdent(theta_z, nam = reTrms),
		pdDiag = pdDiag(theta_z, nam = reTrms),
		pdSymm = pdSymm(theta_z, nam = reTrms),
		pdCompSymm = pdCompSymm(theta_z, nam = reTrms)
)
Psi <- as.matrix(Sigma)/sigma
logDetPsi <- determinant(Psi, logarithm = TRUE)$modulus

Psiinv <- try(solve(Psi), silent = TRUE)
if(inherits(Psiinv, "try-error")) return(Inf)

mmf <- attr(model, "mmf")
mmr <- attr(model, "mmr")

# Multiply design matrices by coefficients
K <- cumsum(sapply(mmf, ncol))
K <- cbind(c(1, K[-length(K)] + 1), K)
for(k in 1:length(mmf)){
	mmf[[k]] <- apply(sweep(mmf[[k]], 2, theta_x[K[k,1]:K[k,2]], "*"), 1, sum)
}
K <- cumsum(sapply(mmr, ncol))
K <- cbind(c(1, K[-length(K)] + 1), K)
for(k in 1:length(mmr)){
	mmr[[k]] <- apply(mmr[[k]]*(z%*%uhat[,K[k,1]:K[k,2]]), 1, sum)
}
# Calculate nonlinear parameters
for(k in 1:length(mmr)) mmf[[names(mmr)[k]]] <- mmr[[k]] + mmf[[names(mmr)[k]]]
# Evaluate model
expr <- deriv3(attr(model, "model"), names(mmf))
eta <- eval(expr, envir = c(data, mmf))
attributes(eta) <- NULL
res <- y - eta

# smoothed loss function
ans <- C_loss_nlqmm(res, uhat, weights, Psiinv, logDetH, logDetPsi, M, N, Q, sigma, tau, omega)
sigma <- ans$hsum/(2*N)
ans <- -ans$val
attr(ans, "resid") <- res
attr(ans, "fitted") <- eta
attr(ans, "sigma") <- sigma
return(ans)
}

#' @export
VarCorr.nlqmm <- function(x, sigma = NULL, ...){

reTrms <- unlist(attr(x$model, "labels")$random)
reTrms <- paste(names(reTrms), reTrms, sep = ".")

Sigma <- switch(x$cov_name,
		pdIdent = pdIdent(x$theta_z, nam = reTrms),
		pdDiag = pdDiag(x$theta_z, nam = reTrms),
		pdSymm = pdSymm(x$theta_z, nam = reTrms),
		pdCompSymm = pdCompSymm(x$theta_z, nam = reTrms)
)

return(Sigma)

}

#' @export
print.nlqmm <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
tau <- x$tau
feTrms <- unlist(attr(x$model, "labels")$fixed)
feTrms <- paste(names(feTrms), feTrms, sep = ".")
reTrms <- unlist(attr(x$model, "labels")$random)
reTrms <- paste(names(reTrms), reTrms, sep = ".")

theta_x <- x$theta_x
sigma <- x$sigma
names(theta_x) <- feTrms

#Sigma <- VarCorr(x)
Sigma <- switch(x$cov_name,
		pdIdent = pdIdent(x$theta_z, nam = reTrms),
		pdDiag = pdDiag(x$theta_z, nam = reTrms),
		pdSymm = pdSymm(x$theta_z, nam = reTrms),
		pdCompSymm = pdCompSymm(x$theta_z, nam = reTrms)
)

cat("Call: ")
dput(x$call)
cat("\n")
cat(paste("Quantile", tau, "\n"))
cat("\n")
cat("Fixed effects:\n")
print.default(format(theta_x, digits = digits), print.gap = 2, 
	quote = FALSE)
cat("\n")
cat("Covariance matrix of the random effects:\n")
print.default(format(as.matrix(Sigma), digits = digits), quote = FALSE)
cat("\n")
cat(paste("Residual scale parameter: ", format(sigma, 
	digits = digits)), "\n")
cat(paste("Log-likelihood:", format(x$logLik, digits = digits), 
	"\n"))
cat(paste("Tuning parameter:", format(x$omega, digits = digits), 
	"\n"))
	
cat(paste("\nNumber of observations:", length(x$y), "\n"))
cat(paste("Number of groups:", x$ngroups, "\n"))

invisible(x)
}

#' @export
ranef.nlqmm <- function(object, ...){

return(object$ranef)

}

#' @export
predict.nlqmm <- function(object, level = 0, newdata = NULL, ...){

mmf <- attr(object$model, "mmf")
theta <- object$theta_x

K <- cumsum(sapply(mmf, ncol))
K <- cbind(c(1, K[-length(K)] + 1), K)
for(k in 1:length(mmf)){
	mmf[[k]] <- apply(sweep(mmf[[k]], 2, theta[K[k,1]:K[k,2]], "*"), 1, sum)
	mmf[[k]] <- mmf[[k]][object$revOrder]
}

model <- attr(object$model, "model")[c(1,3)]

val <- eval(deriv(model, "x"), envir = c(object$mfArgs$data,mmf))
attr(val, "gradient") <- NULL
return(val)

}

bootData <- function(x, group, numeric = TRUE){

group.old <- group
	if(is.factor(group.old)){
		if(numeric){
			group <- as.numeric(levels(group))[group]
		} else {
			group <- as.character(levels(group))[group]
		}
	}
x <- x[order(group),]
group <- sort(group)
grp <- unique(group)
M <- length(grp)
n <- table(group)

id <- sort(sample(grp, M, replace = TRUE))
dd <- do.call("rbind", lapply(id, function(i) subset(x, group == i)))

nid <- rep(grp, n[id])
if(is.factor(group.old)){
	nid <- factor(nid, levels(group.old))
}
dd$newid <- nid
dd
}
############################################################################
