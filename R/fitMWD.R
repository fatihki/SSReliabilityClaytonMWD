#' Classical estimations for the parameters of the Modified Weibull Distribution (MWD) 
#' @title Estimating parameters of the Modified Weibull Distribution (MWD) 
#' 
#' @description
#' Estimates the parameters of the Modified Weibull Distribution (MWD)
#' using classical methods.
#'
#' @import stats
#' 
#' @name fitMWD
#'  
#' @param data Vector of observations.
#' @param est.method Used method for estimating the parameters such as maximum likelihood estimate "MLE",
#' the least square estimation "LSE", the weighted least square estimation "WLSE", and
#' the maximum product of spacing estimates "MPS".
#' @param opt.method The optimization method for \code{optim} function such as "Nelder-Mead", "BFGS", "CG", 
#' "L-BFGS-B", "SANN" and "Brent" to be used for estimating the parameters.
#' @param starts Initial values for the parameters to be optimized over.
#' @param lower Numeric vector specifying the lower bounds of the parameters
#' for bounded optimization methods (e.g., \code{"L-BFGS-B"}). Ignored if \code{NULL}.
#' @param upper Numeric vector specifying the upper bounds of the parameters
#' for bounded optimization methods. Ignored if \code{NULL}.
#' @param verbose Logical; if \code{TRUE}, progress and intermediate
#' results from the optimization procedure are printed. Default is \code{FALSE}.
#' @param ... Further arguments to be passed to \code{optim} function such as lower and upper limits, hessian, etc.
#' 
#' 
#' @details 
#' The Modified Weibull Distribution (Lai et al., 2003) has cumulative distribution
#' function (CDF) and probability density function (PDF) given by
#'
#' \deqn{
#' F(x) = 1 - \exp\left(-a x^b \exp(\lambda x)\right),
#' }
#' \deqn{
#' f(x) = a (b + \lambda x) x^{b - 1} \exp(\lambda x)
#' \exp\left(-a x^b \exp(\lambda x)\right),
#' }
#'
#' where \eqn{x > 0}, \eqn{a > 0} is a scale parameter, \eqn{b \ge 0} is a shape parameter,
#' and \eqn{\lambda \ge 0} is a flexibility parameter controlling the growth rate of the hazard function.
#' 
#' The model parameters are estimated using several classical methods:
#'
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):}
#'   Obtains parameter estimates by maximizing the log-likelihood
#'   function of the observed data under the assumed MWD distribution.
#'
#'   \item \strong{Least Squares Estimation (LSE):}
#'   Estimates parameters by minimizing the sum of squared differences
#'   between the empirical distribution function and theoretical distribution function
#'   of MWD. 
#'   The Benard's approximation is used for the empirical distribution function
#'    at the ordered observations as \eqn{F(x_(i)) = (i-0.3)/(n+0.4), i=1,...,n.}.
#'
#'   \item \strong{Weighted Least Squares Estimation (WLSE):}
#'   A modification of LSE that assigns weights to the squared differences.
#'   The weights are taken as \eqn{w_i = ((n+1)^2(n+2)) / (i(n-i+1)), i=1,...,n. }
#'
#'   \item \strong{Maximum Product of Spacings (MPS):}
#'   Estimates parameters by maximizing the product of spacings between
#'   consecutive values of the fitted distribution function, providing
#'   a robust alternative to MLE, particularly in small samples.
#' }
#' Further details can be found in Kizilaslan (2026).
#'
#' @return A list containing:
#' \item{estimates}{Estimated values of the model parameters.}
#' \item{measures}{Model selection criteria, including the log-likelihood, AIC, and BIC, 
#' evaluated at the estimated parameter values.}
#' \item{initials}{Initial values used in the optimization procedure.}
#' \item{opt.fit}{Full optimization output.}
#' 
#' @references
#' Lai, C. D., Xie, M., and Murthy, D. N. P. (2003).
#' \href{A modified Weibull distribution.}{https://doi.org/10.1109/TR.2002.805788}
#' \emph{IEEE Transactions on Reliability}, \strong{52}(1), 33--37.
#'
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' arXiv preprint. Available at
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}.
#' 
#' @examples 
#' # generate data from MWD(a, b, lambda)
#' n <- 100
#' a <- 0.75; b <- 1.25; lambda <- 0.60
#' set.seed(123)
#' dat <- rMweibull(n, a, b, lambda)
#' init <- runif(3)
#' 
#' # Fit MWD to dat.
#' fit.mle <- fitMWD(data = dat, est.method = "mle", opt.method = "L-BFGS-B", starts = init,
#'                   lower = c(1e-05,1e-05,1e-05), upper = c(Inf,Inf,Inf), hessian = FALSE )
#' fit.mle$estimates
#' 
#' fit.lse <- fitMWD(data = dat, est.method = "lse", opt.method = "L-BFGS-B", starts = init,
#'                   lower = c(1e-05,1e-05,1e-05), upper = c(Inf,Inf,Inf), hessian = FALSE )
#' fit.lse$estimates
#'
#' fit.wlse <- fitMWD(data = dat, est.method = "wlse", opt.method = "L-BFGS-B", starts = init,
#'                    lower = c(1e-05,1e-05,1e-05), upper = c(Inf,Inf,Inf), hessian = FALSE )
#' fit.wlse$estimates
#' 
#' fit.mps <- fitMWD(data = dat, est.method = "mps", opt.method = "L-BFGS-B", starts = init,
#'                   lower = c(1e-05,1e-05,1e-05), upper = c(Inf,Inf,Inf), hessian = FALSE )
#' fit.mps$estimates
#'   
#' @export
fitMWD <- function(data, est.method, opt.method, starts, lower = NULL, upper = NULL, verbose = FALSE, ... ){

  if (!est.method %in% c("mle", "lse", "wlse", "mps"))
    stop("Estimation method class misspelled. Please check it.")
  
  if (!opt.method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
    stop("Optimization method class misspelled. Please check it.")
  

  if(est.method == "mle"){
    
    opt.args <- list(
      par     = starts,
      fn      = nll_MWD,
      gr      = grad_nll_MWD,
      x       = data,
      method  = opt.method,
      control = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
      out <- tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                         ),
    error = function(e) {
      if (verbose) {
      message("Optimization failed") }
      NULL
      }
      )
      
    }
  
  if(est.method == "lse"){
    
    opt.args <- list(
      par     = starts,
      fn      = lse_MWD,
      gr      = grad_lse_MWD,
      x       = data,
      method  = opt.method,
      control = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
    out <-  tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                        ),
    error = function(e) {
      if (verbose) {
      message("Optimization failed")
      }
      NULL
    }
    )

  }
  
  if(est.method == "wlse"){
    
    opt.args <- list(
      par     = starts,
      fn      = wlse_MWD,
      gr      = grad_wlse_MWD,
      x       = data,
      method  = opt.method,
      control = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
    out <-  tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                        ),
    error = function(e) {
      if (verbose) {
      message("Optimization failed")
      }
      NULL
    }
    )

  }
  
  if(est.method == "mps"){
    
    opt.args <- list(
      par     = starts,
      fn      = mps_MWD,
      gr      = grad_mps_MWD,
      x       = data,
      method  = opt.method,
      control = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
    out <-  tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                        ),
    error = function(e) {
      if (verbose) {
      message("Optimization failed")
      }
      NULL
    }
    )

  }
  
  if (is.null(out)) {
    if(verbose){
    message("Optimization failed -- exiting this run.")
    }
    return(NULL)  # or stop() if you want to terminate entirely
  }
  
  par.est <- out$par
  names(par.est) <- c("a","b","lambda")
  n <- length(data)
  log.likelihod <- -1*nll_MWD(par.est, data)
  AIC <- -2*log.likelihod + 2*length(par.est)
  BIC <- -2*log.likelihod + log(n)*length(par.est)
  out.measures <- cbind(log.likelihod, AIC, BIC)
  colnames(out.measures) <- c("log.likelihood", "AIC","BIC")
  
  return( list("estimates" = par.est, "measures"= out.measures, initials = starts, opt.fit = out) )
}
#' 
#' @keywords internal
nll_MWD <- function(par, x){
  a <- par[1]; b <- par[2]; lambda <- par[3]
  if(a <= 0 || b < 0 || lambda < 0 ) return(-Inf)  
  logf <- log(a) + (b-1)*log(x) + lambda*x + log(b + lambda*x) - a*x^b*exp(lambda*x)
  return(-sum(logf))
}
#' @keywords internal
grad_nll_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  if(a <= 0 || b < 0  || lambda < 0 ) return(rep(NA,3))  # positivity
  n <- length(x)
  xb <- x^b
  exp_lx <- exp(lambda*x)
  
  grad_a <- -sum(1/a - xb*exp_lx)
  grad_b <- -sum(log(x) + 1/(b + lambda*x) - a * xb * exp_lx * log(x))
  grad_lambda <- -sum(x + x/(b + lambda*x) - a * xb * exp_lx * x)
  
  return( c(grad_a, grad_b, grad_lambda) )
}
#'
#' @keywords internal
lse_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} ) 
  Fmodel <- 1 - exp(-a * x^b * exp(lambda*x))
  return( sum( (Fmodel - Fhat)^2) )
}
#'
#' @keywords internal
grad_lse_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} )
  Fmodel <- 1 - exp(-a * x^b * exp(lambda*x))
  grad_Fa <- (1-Fmodel) * x^b * exp(lambda*x)
  grad_Fb <- (1-Fmodel) * a * x^b * log(x) * exp(lambda*x)
  grad_Flambda <- (1-Fmodel) * a * x^(b+1) * exp(lambda*x)

  grad_a <- sum( 2*(Fmodel - Fhat)*grad_Fa )
  grad_b <- sum( 2*(Fmodel - Fhat)*grad_Fb )
  grad_lambda <- sum( 2*(Fmodel - Fhat)*grad_Flambda )

  return( c(grad_a, grad_b, grad_lambda) )
}
#' 
#' @keywords internal
wlse_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} )
  Fmodel <- 1 - exp(-a * x^b * exp(lambda*x))
  w     <- sapply( 1:n, function(i){ ( (n+1)^2*(n+2) ) / ( i*(n-i+1) )} )
  return( sum( w*(Fmodel - Fhat)^2) )
}
#'
#' @keywords internal
grad_wlse_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} ) 
  Fmodel <- 1 - exp(-a * x^b * exp(lambda*x))
  w     <- sapply( 1:n, function(i){ ( (n+1)^2*(n+2) ) / ( i*(n-i+1) )} )

  grad_Fa <- (1-Fmodel) * x^b * exp(lambda*x)
  grad_Fb <- (1-Fmodel) * a * x^b * log(x) * exp(lambda*x)
  grad_Flambda <- (1-Fmodel) * a * x^(b+1) * exp(lambda*x)

  grad_a <- sum( 2*w*(Fmodel - Fhat)*grad_Fa )
  grad_b <- sum( 2*w*(Fmodel - Fhat)*grad_Fb )
  grad_lambda <- sum( 2*w*(Fmodel - Fhat)*grad_Flambda )

  return( c(grad_a, grad_b, grad_lambda) )
}
#' 
#' @keywords internal
mps_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fmodel <- 1-exp(-a * x^b * exp(lambda*x))
  D     <- diff( c(0, Fmodel, 1) )
  return( -sum(log(D))/(n+1) )
}
#'
#' @keywords internal
grad_mps_MWD <- function(par, x) {
  a <- par[1]; b <- par[2]; lambda <- par[3]
  x    <- sort(x)
  n    <- length(x)
  Fmodel <- 1 - exp(-a * x^b * exp(lambda*x))
  D      <- diff( c(0, Fmodel, 1) )

  grad_Fa <- (1-Fmodel) * x^b * exp(lambda*x)
  D_grad_Fa <- diff( c(0, grad_Fa, 0) )
  grad_Fb <- (1-Fmodel) * a * x^b * log(x) * exp(lambda*x)
  D_grad_Fb <- diff( c(0, grad_Fb, 0) )
  grad_Flambda <- (1-Fmodel) * a * x^(b+1) * exp(lambda*x)
  D_grad_Flambda <- diff( c(0, grad_Flambda, 0) )

  grad_a <- -sum( (D_grad_Fa/D)/(n+1)  )
  grad_b <- -sum( (D_grad_Fb/D)/(n+1)   )
  grad_lambda <- -sum( (D_grad_Flambda/D)/(n+1)  )

  return( c(grad_a, grad_b, grad_lambda) )
}
