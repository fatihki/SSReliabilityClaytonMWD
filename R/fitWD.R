#' Classical estimations fo the parameters of the two-parameter Weibull distribution
#' @title Estimating parameters of the two-parameter Weibull distribution
#' 
#' @description
#' Estimates the parameters of the two-parameter Weibull distribution
#' using classical methods.
#' 
#' @name fitWD
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
#' The two-parameter Weibull Distribution has cumulative distribution
#' function (CDF) and probability density function (PDF) given by
#'
#' \deqn{F(x) = 1 - \exp(-a x^b),}
#' \deqn{f(x) = a b x^{b - 1} \exp(-a x^b ),}
#' 
#' where \eqn{x > 0}, \eqn{a > 0} is the scale parameter and \eqn{b > 0} is the shape parameter.
#'
#' The model parameters are estimated using several classical methods:
#'
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):}
#'   Obtains parameter estimates by maximizing the log-likelihood
#'   function of the observed data under the assumed two-parameter Weibull distribution.
#'
#'   \item \strong{Least Squares Estimation (LSE):}
#'   Estimates parameters by minimizing the sum of squared differences
#'   between the empirical distribution function and theoretical distribution function
#'   of MWD. 
#'   The Benard’s approximation is used for the empirical distribution function
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
#' 
#' @return A list containing:
#' \item{estimates}{Estimated values of the model parameters.}
#' \item{measures}{Model selection criteria, including the log-likelihood, AIC, and BIC, 
#' evaluated at the estimated parameter values.}
#' \item{initials}{Initial values used in the optimization procedure.}
#' \item{opt.fit}{Full optimization output.}
#' 
#' @examples 
#' # generate data from WD(a, b)
#' n <- 50
#' a <- 0.75; b <- 1.25; lambda <- 0 # reduces two-parameter Weibull distribution
#' set.seed(123)
#' X <- rMweibull(n, a, b, lambda)
#' init <- runif(2)
#' 
#' # fit model
#' fit.mle <- fitWD(data = X, est.method = "mle", opt.method = "L-BFGS-B", starts = init,
#'                  lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
#' fit.mle$estimates
#' 
#' fit.lse <- fitWD(data = X, est.method = "lse", opt.method = "L-BFGS-B", starts = init,
#'                  lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
#' fit.lse$estimates
#' 
#' fit.wlse <- fitWD(data = X, est.method = "wlse", opt.method = "L-BFGS-B", starts = init,
#'                   lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
#' fit.wlse$estimates
#' 
#' fit.mps <- fitWD(data = X, est.method = "mps", opt.method = "L-BFGS-B", starts = init,
#'                  lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
#' fit.mps$estimates
#' 
#' @export
fitWD <- function(data, est.method, opt.method, starts, lower = NULL, upper = NULL, verbose = FALSE, ... ){
  
  if (!est.method %in% c("mle", "lse", "wlse", "mps"))
    stop("Estimation method class misspelled. Please check it.")
  
  if (!opt.method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "optimx"))
    stop("Optimization method class misspelled. Please check it.")
  
  
  if(est.method == "mle"){
    
    opt.args <- list(
      par     = starts,
      fn      = nll_WD,
      gr      = grad_nll_WD,
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
      message("Optimization failed")
      }
      NULL
    }
    )
  }
  
  if(est.method == "lse"){
    
    opt.args <- list(
      par     = starts,
      fn      = lse_WD,
      gr      = grad_lse_WD,
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
      fn      = wlse_WD,
      gr      = grad_wlse_WD,
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
      fn      = mps_WD,
      gr      = grad_mps_WD,
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
    if (verbose) {
    message("Optimization failed — exiting this run.")
    }
    return(NULL)  # or stop() if you want to terminate entirely
  }
  
  
  par.est <- out$par
  names(par.est) <- c("a","b")
  n <- length(data)
  log.likelihod <- -1*nll_WD(par.est, data)
  AIC <- -2*log.likelihod + 2*length(par.est)
  BIC <- -2*log.likelihod + log(n)*length(par.est)
  out.measures <- cbind(log.likelihod, AIC, BIC)
  colnames(out.measures) <- c("log.likelihood", "AIC","BIC")
  
  return( list("estimates" = par.est, "measures"= out.measures, initials = starts, opt.fit = out) )
}
#'
#' @keywords internal
nll_WD <- function(par, x){
  a <- par[1]; b <- par[2];
  if(a <= 0 || b < 0 ) return(-Inf)  
  logf <- log(a) + (b-1)*log(x) + log(b) - a*x^b
  return(-sum(logf))
}
#'
#' @keywords internal
grad_nll_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  if(a <= 0 || b < 0   ) return(rep(NA,2))  
  n <- length(x)
  xb <- x^b
  
  grad_a <- -sum(1/a - xb)
  grad_b <- -sum(log(x) + 1/(b ) - a * xb * log(x))
  
  return( c(grad_a, grad_b) )
}
#' 
#' @keywords internal
lse_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} ) 
  Fmodel <- 1 - exp(-a * x^b  )
  return( sum( (Fmodel - Fhat)^2) )
}
#'
#' @keywords internal
grad_lse_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} ) 
  Fmodel <- 1 - exp(-a * x^b)
  grad_Fa <- (1-Fmodel) * x^b 
  grad_Fb <- (1-Fmodel) * a * x^b * log(x) 
  
  grad_a <- sum( 2*(Fmodel - Fhat)*grad_Fa )
  grad_b <- sum( 2*(Fmodel - Fhat)*grad_Fb )
  
  return( c(grad_a, grad_b) )
}
#'
#' @keywords internal
wlse_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} )
  Fmodel <- 1 - exp(-a * x^b )
  w     <- sapply( 1:n, function(i){ ( (n+1)^2*(n+2) ) / ( i*(n-i+1) )} )
  return( sum( w*(Fmodel - Fhat)^2) )
}
#' 
#' @keywords internal
grad_wlse_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fhat <- sapply(1:n, function(i){ (i-0.3) / (n+0.4)} )
  Fmodel <- 1 - exp(-a * x^b)
  w     <- sapply( 1:n, function(i){ ( (n+1)^2*(n+2) ) / ( i*(n-i+1) )} )
  
  grad_Fa <- (1-Fmodel) * x^b 
  grad_Fb <- (1-Fmodel) * a * x^b * log(x) 
  
  grad_a <- sum( 2*w*(Fmodel - Fhat)*grad_Fa )
  grad_b <- sum( 2*w*(Fmodel - Fhat)*grad_Fb )
  
  return( c(grad_a, grad_b) )
}
#' 
#' @keywords internal
mps_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fmodel <- 1-exp(-a * x^b)
  D     <- diff( c(0, Fmodel, 1) )
  return( -sum(log(D))/(n+1) )
}
#'
#' @keywords internal
grad_mps_WD <- function(par, x) {
  a <- par[1]; b <- par[2]
  x    <- sort(x)
  n    <- length(x)
  Fmodel <- 1 - exp(-a * x^b )
  D      <- diff( c(0, Fmodel, 1) )
  
  grad_Fa <- (1-Fmodel) * x^b 
  D_grad_Fa <- diff( c(0, grad_Fa, 0) )
  grad_Fb <- (1-Fmodel) * a * x^b * log(x)
  D_grad_Fb <- diff( c(0, grad_Fb, 0) )
  
  grad_a <- -sum( (D_grad_Fa/D)/(n+1)  )
  grad_b <- -sum( (D_grad_Fb/D)/(n+1)   )
  
  return( c(grad_a, grad_b) )
}
