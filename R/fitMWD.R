#' Fit the Modified Weibull Distribution (MWD)
#'
#' @title Estimation of Parameters for the Modified Weibull Distribution
#' 
#' @description
#' Estimates the parameters of the Modified Weibull Distribution (MWD)
#' using classical estimation methods.
#'
#' @import stats
#' 
#' @name fitMWD
#'  
#' @param data Numeric vector of observations.
#'
#' @param est.method Character string specifying the estimation method.
#' Options include \code{"MLE"}, \code{"LSE"}, \code{"WLSE"}, and \code{"MPS"}.
#'
#' @param opt.method Character string specifying the optimization method
#' used in \code{optim}, such as \code{"Nelder-Mead"}, \code{"BFGS"},
#' \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}, or \code{"Brent"}.
#'
#' @param starts Numeric vector of initial values for the parameters
#' 
#' @param lower Numeric vector of lower bounds for parameters in constrained optimization.
#' Ignored if \code{NULL}.
#'
#' @param upper Numeric vector of upper bounds for parameters in constrained optimization.
#'
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#'
#' @param ... Additional arguments passed to \code{optim}.
#' 
#' 
#' @details
#' The Modified Weibull Distribution (Lai et al., 2003) has cumulative
#' distribution function (CDF) and probability density function (PDF):
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
#' The parameters are estimated using the following methods:
#'
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):}
#'   Maximizes the log-likelihood under the MWD model.
#'
#'   \item \strong{Least Squares Estimation (LSE):}
#'   Minimizes squared differences between empirical and theoretical CDFs.
#'   The empirical CDF uses Benard's approximation:
#'   \eqn{F(x_{(i)}) = (i - 0.3)/(n + 0.4)}, for \eqn{i = 1, \dots, n}.
#'
#'   \item \strong{Weighted Least Squares Estimation (WLSE):}
#'   A modification of LSE that assigns weights to the squared differences.
#'   Uses weights
#'   \eqn{w_i = \frac{(n+1)^2(n+2)}{i(n-i+1)}}, for \eqn{i = 1, \dots, n}.
#'   
#'   \item \strong{Maximum Product of Spacings (MPS):}
#'   Maximizes the product of spacings of the fitted CDF.
#' }
#' 
#' Further details are provided in Kizilaslan (2026).
#'
#' @return A list containing:
#' \item{estimates}{Named numeric vector of estimated parameters \eqn{(a, b, \lambda)}.}
#' \item{measures}{Numeric vector of model selection criteria (log-likelihood, AIC, BIC).}
#' \item{initials}{Initial values used in the optimization.}
#' \item{opt.fit}{Full output from \code{optim}.}
#' 
#' @references
#' Lai, C. D., Xie, M., and Murthy, D. N. P. (2003).
#' \href{https://doi.org/10.1109/TR.2002.805788}{A modified Weibull distribution.}
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
