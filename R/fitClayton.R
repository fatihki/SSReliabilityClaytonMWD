#' Estimate the Clayton Copula Parameter
#'
#' @title Estimation of the Clayton Copula Dependence Parameter
#'
#' @description
#' Estimates the dependence parameter \eqn{\theta} of the Clayton copula
#' based on observed data from a stress–strength model.
#' 
#' @import stats
#'  
#' @name fitClayton
#'
#' @param x Numeric vector. Observations of the strength variable \eqn{X}.
#' @param y Numeric vector. Observations of the stress variable \eqn{Y}.
#' 
#' @param est.method Character string specifying the estimation method used.
#'  Options include \code{"MLE"}, \code{"LSE"}, \code{"WLSE"}, and \code{"MPS"}.
#'  
#' @param opt.method Character string specifying the optimization method used in \code{optim}.
#' Common options include \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"},
#' \code{"L-BFGS-B"}, \code{"SANN"}, and \code{"Brent"}.
#' 
#' @param start Numeric scalar. Initial value for \eqn{\theta}.
#' 
#' @param estimates A named list of estimated marginal parameters:
#' \eqn{(a_1, b_1, \lambda_1)} for strength and
#' \eqn{(a_2, b_2, \lambda_2)} for stress.
#' 
#' @param lower Numeric vector. Lower bounds for parameters in constrained optimization.
#' Only used if supported by \code{opt.method}.
#'
#' @param upper Numeric vector. Upper bounds for parameters in constrained optimization.
#' Only used if supported by \code{opt.method}.
#' 
#' @param verbose Logical; if \code{TRUE}, progress and intermediate
#' results from the optimization procedure are printed. Default is \code{FALSE}.
#'
#' @param ... Additional arguments passed to \code{optim}.
#' 
#' @details 
#' The Clayton copula is defined as
#' \deqn{
#' C(u,v;\theta) = \left(u^{-\theta} + v^{-\theta} - 1\right)^{-1/\theta},
#' }
#' where \eqn{\theta > 0}.
#'
#' The parameter is estimated using the following methods:
#'
#' \itemize{
#'   \item \strong{Maximum Likelihood Estimation (MLE):}
#'   Maximizes the joint log-likelihood under the assumed model.
#'
#'   \item \strong{Least Squares Estimation (LSE):}
#'   Minimizes squared differences between empirical and theoretical CDFs.
#'   The empirical CDF uses Benard's approximation:
#'   \eqn{F(x_{(i)}) = (i - 0.3)/(n + 0.4)}, for \eqn{i = 1, \dots, n}.
#'
#'   \item \strong{Weighted Least Squares Estimation (WLSE):}
#'   Uses weights
#'   \eqn{w_i = \frac{(n+1)^2(n+2)}{i(n-i+1)}}, for \eqn{i = 1, \dots, n}.
#'
#'   \item \strong{Maximum Product of Spacings (MPS):}
#'   Maximizes the product of spacings of the fitted distribution function,
#'   providing a robust alternative to MLE.
#' }
#'
#' Further theoretical details are provided in Kizilaslan (2026).
#' 
#' @return A list containing:
#' \item{estimate}{Estimate of the Clayton copula parameter, \eqn{\theta}.}
#' \item{opt.fit}{Full optimization result.}
#' 
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' arXiv preprint. Available at
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}.
#' 
#' @export
fitClayton <- function(x, y, est.method, opt.method, start, estimates, 
                       lower = NULL, upper = NULL, verbose = FALSE, ... ){
  
  if (!est.method %in% c("mle", "lse", "wlse", "mps"))
    stop("Estimation method class misspelled. Please check it.")
  
  if (!opt.method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
    stop("Optimization method class misspelled. Please check it.")
  
  if(est.method == "mle"){
    
    opt.args    <- list(
      par       = start,
      fn        = nll_clayton,
      gr        = grad_nll_clayton,
      x         = x,
      y         = y, 
      estimates = estimates,
      method    = opt.method,
      control   = list(maxit = 5000)
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
    
    opt.args    <- list(
      par       = start,
      fn        = lse_clayton,
      x         = x,
      y         = y, 
      estimates = estimates,
      method    = opt.method,
      control   = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
      out <-  tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                          ),
      error = function(e) {
        if (verbose) {
          message("Optimization failed") }
        NULL
      }
      )

  }
  
  if(est.method == "wlse"){
    
    opt.args    <- list(
      par       = start,
      fn        = wlse_clayton,
      x         = x,
      y         = y, 
      estimates = estimates,
      method    = opt.method,
      control   = list(maxit = 5000)
    )
    
    if (opt.method == "L-BFGS-B") {
      opt.args$lower <- lower
      opt.args$upper <- upper
    }
    
      out <-  tryCatch( suppressWarnings( do.call(optim, c(opt.args, list(...)))
                                          ),
      error = function(e) {
        if (verbose) {
          message("Optimization failed") }
        NULL
      }
      )

  }
  
  if(est.method == "mps"){
    
    opt.args    <- list(
      par       = start,
      fn        = mps_clayton,
      x         = x,
      y         = y, 
      estimates = estimates,
      method    = opt.method,
      control   = list(maxit = 5000)
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
  
  if (is.null(out)) {
    if(verbose){
      message("Optimization failed -- exiting this run.")
    }
    return(NULL)  # or stop() if you want to terminate entirely
  }
  
  par.est <- out$par
  n <- length(x)
  log.likelihod <- -nll_clayton(out$par, x, y, estimates)
  AIC <- -2*log.likelihod + 2*length(out$par)
  BIC <- -2*log.likelihod + log(n)*length(out$par)
  out.measures <- cbind(log.likelihod, AIC, BIC)
  colnames(out.measures) <- c("log.likelihood", "AIC","BIC")
  
  return( list("estimates" = par.est, "measures"= out.measures, opt.fit = out) )
}
#'
# ---------------------------------------------------------------------------
# Objective functions are used in optimization of the methods in fitClayton
# ---------------------------------------------------------------------------
#' @keywords internal
nll_clayton  <- function(par, x, y, estimates){ 
  if(par < 0 ) return(-Inf)  
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  logL <- n*log(1+par) - (1+par)*sum(log(u*v)) - (2+(1/par))*sum( log(u^(-par) + v^(-par) -1) )
  return(-logL)
}
#' @keywords internal 
grad_nll_clayton  <- function(par, x, y, estimates){ 
  if(par < 0 ) return(-Inf)  
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  grad_par1 <- n/(1+par) - sum(log(u*v)) + (1/par^2)*sum( log(u^(-par) + v^(-par) -1) ) 
  grad_par2 <- (2+(1/par))*sum( ( u^(-par)*log(u) + v^(-par)*log(v) )/(u^(-par) + v^(-par) -1) )
  
  return(-(grad_par1+grad_par2) )
}
#'
#' the derivative of l3(theta) wrt theta as above "grad_nll_clayton" used in the next function "grad2_nll3"
#'@noRd 
grad_nll3  <- function(par, x, y){ 
  if(any(par< 0)) return(-Inf)  
  n <- length(x)
  a1 <- par[1]; b1 <- par[2]; lambda1 <- par[3]; 
  a2 <- par[4]; b2 <- par[5]; lambda2 <- par[6]; theta <- par[7]
  u <- pMweibull(x,a1,b1,lambda1)
  v <- pMweibull(y,a2,b2,lambda2)

  grad_par1 <- n/(1+theta) - sum(log(u*v)) + (1/theta^2)*sum( log(u^(-theta) + v^(-theta) -1) ) 
  grad_par2 <- (2+(1/theta))*sum( ( u^(-theta)*log(u) + v^(-theta)*log(v) )/(u^(-theta) + v^(-theta) -1) )
  return(-(grad_par1+grad_par2) )
}
#'
#' 2nd order derivatives of grad_nll3 wrt parameter a1,b1,lambda1,a2,b2,lambda2,theta respectively.
#' It is used for ACI of the parameters based on MLEs.
#' @noRd
grad2_nll3 <- function(par, x, y){ 
  a1 <- par[1]; b1 <- par[2]; lambda1 <- par[3]; 
  a2 <- par[4]; b2 <- par[5]; lambda2 <- par[6]; theta <- par[7]
     grad_nll3_r1 <- function(par){
          return(grad_nll3(par, x, y)) 
       }
  grad2_l3 <- numDeriv::grad( grad_nll3_r1, x = par )
  return(grad2_l3)
}
#'
#' @noRd 
empirical_cdf <- function(u, v){
  n <- length(u)
  H <- c()
  for(i in 1:n){
    H[i] <- mean(u <= u[i] & v <= v[i])
    }
  return(H)
}
#'
#' @keywords internal 
lse_clayton <- function(par, x, y, estimates) {
  
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Hhat <- empirical_cdf(u, v)
  Cval <- Clayton_Copula(u, v, par)
  
  return( sum((Cval - Hhat)^2) )
}
#'
#' @keywords internal
wlse_clayton <- function(par, x, y, estimates) {
  
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Hhat <- empirical_cdf(u, v)
  w <- 1 / (Hhat * (1 - Hhat) +  1e-6)   # Variance-stabilizing weights
  Cval <- Clayton_Copula(u, v, par)
  
  return( sum( w* (Cval - Hhat)^2 ) )
}
#' 
#' @keywords internal
mps_clayton <- function(par, x, y, estimates) {
  if(par < 0 ) return(-Inf)  
  n <- length(x) 
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Ccopula <- Clayton_Copula(u, v, par)
  Ccopula <- sort(Ccopula)
  D3 <- diff( c(0, Ccopula, 1) )
  # multiplying with -1 for minimizing
  return( -sum(log(D3))/(n+1) )
}
#'
#' Kendall's Tau Estimator for the Clayton Copula Parameter
#'
#' @title Kendall's Tau-based Estimation of the Clayton Copula Parameter
#' 
#' @description
#' Estimates the dependence parameter \eqn{\theta} of the Clayton copula
#' using Kendall's tau-based moment estimator.
#' 
#' @import stats
#'  
#' @name theta_Ktau_estimate
#' 
#' @param data A list containing two numeric vectors:
#'   \code{X} (strength) and \code{Y} (stress).
#'
#' @return A numeric scalar giving the estimate of \eqn{\theta}
#' based on Kendall's tau (\eqn{\tau}).
#'
#' @details
#' The estimator is derived from the relationship between Kendall's tau
#' and the Clayton copula parameter:
#' \eqn{\tau = \theta / (\theta + 2)}.
#' 
#' @examples
#' set.seed(123)
#' n <- 50
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 5 # 1, 2, 3, 4
#' # data generation
#' dat <- SSReliabilityClaytonMWD::rMweibull_Clayton(n, a1, b1, lambda1, a2, b2, lambda2, theta)
#' theta_Ktau_estimate(dat)
#' 
#' @export
theta_Ktau_estimate <-function(data){
  tau_hat <- cor(data$X, data$Y, method = "kendall")
  theta_tau <- (2 * tau_hat) / (1 - tau_hat)
  return(theta_tau)
}
#'
# -------------------------------
# One-step LSE estimate
# -------------------------------
#' 
#' One-Step LSE Estimator for the Clayton Copula Parameter
#'
#' @title One-Step Least Squares Estimation of the Clayton Copula Parameter
#' 
#' @description
#' Computes a one-step least squares estimator (LSE) of the Clayton copula
#' dependence parameter \eqn{\theta}. The estimator is obtained via a
#' second-order Taylor expansion of the Clayton copula \eqn{C_{\theta}(u, v)}
#' around an initial value \eqn{\theta_0}, typically the Kendall's
#' tau-based moment estimate.
#'
#' @name LSE_clayton_onestep
#' 
#' @param par Numeric scalar. Initial estimate of \eqn{\theta}, typically
#' obtained from Kendall's tau.
#'
#' @param x Numeric vector. Observations of the strength variable \eqn{X}.
#'
#' @param y Numeric vector. Observations of the stress variable \eqn{Y}.
#'
#' @param estimates A named list of marginal parameter estimates:
#' \eqn{(a_1, b_1, \lambda_1)} for strength and
#' \eqn{(a_2, b_2, \lambda_2)} for stress.
#' 
#' @details
#' The one-step estimator is constructed by substituting a second-order Taylor
#' expansion of the Clayton copula \eqn{C_{\theta}(u, v)} into the least
#' squares estimating equation, evaluated at \eqn{\theta_0}, and solving 
#' analytically for \eqn{\theta}.This avoids iterative numerical optimisation 
#' and yields a closed-form estimation of \eqn{\theta}.
#'
#' Further theoretical details are provided in Kizilaslan (2026).
#' 
#' @return Numeric scalar. One-step LSE estimate of the dependence parameter \eqn{\theta}.
#' 
#' @seealso
#' \code{\link{WLSE_clayton_onestep}} for the weighted LSE version,
#' \code{\link{theta_Ktau_estimate}} for the Kendall's Tau-based estimate.
#'
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @export
LSE_clayton_onestep <- function(par, x, y, estimates) {
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Hhat <- empirical_cdf(u, v)
  
  # Compute derivatives at par = theta
  C_vals <- Clayton_Copula(u, v, par)
  dC_vals <- dClayton_Copula_theta(u, v, par)
  d2C_vals <- d2Clayton_Copula_theta(u, v, par)
  
  # One-step formula for LSE
  B <- 2 * sum( (C_vals - Hhat) * dC_vals )
  C_val <- sum( dC_vals^2 + (C_vals - Hhat) * d2C_vals )
  
  theta_hat <- par - B / (2*C_val)
  
  return(theta_hat)
}
# -------------------------------
# One-step WLSE estimate
# -------------------------------
#' One-Step WLSE Estimator for the Clayton Copula Parameter
#'
#' @title One-Step Weighted Least Squares Estimation of the Clayton Copula Parameter
#' 
#' @description
#' Computes a one-step weighted least squares estimator (WLSE) of the Clayton copula
#' dependence parameter \eqn{\theta}. The estimator is obtained via a
#' second-order Taylor expansion of the Clayton copula \eqn{C_{\theta}(u, v)}
#' around an initial value \eqn{\theta_0}, typically the Kendall's
#' tau-based moment estimate.
#' 
#' @name WLSE_clayton_onestep
#' 
#' @param par Numeric scalar. Initial estimate of \eqn{\theta}, typically
#' obtained from Kendall's tau.
#'
#' @param x Numeric vector. Observations of the strength variable \eqn{X}.
#'
#' @param y Numeric vector. Observations of the stress variable \eqn{Y}.
#'
#' @param estimates A named list of marginal parameter estimates:
#' \eqn{(a_1, b_1, \lambda_1)} for strength and
#' \eqn{(a_2, b_2, \lambda_2)} for stress.
#' 
#' @details
#' The one-step estimator is constructed by substituting a second-order Taylor
#' expansion of the Clayton copula \eqn{C_{\theta}(u, v)} into the weighted least
#' squares estimating equation, evaluated at \eqn{\theta_0}, and solving
#' analytically for \eqn{\theta}.This avoids iterative numerical optimisation 
#' and yields a closed-form estimation of \eqn{\theta}.
#'
#' Further theoretical details are provided in Kizilaslan (2026).
#' 
#' @return Numeric scalar. One-step WLSE estimate of the dependence parameter \eqn{\theta}.
#' 
#' @seealso
#' \code{\link{LSE_clayton_onestep}} for the LSE version,
#' \code{\link{theta_Ktau_estimate}} for the Kendall's Tau-based estimate.
#'
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @export
WLSE_clayton_onestep <- function( par, x, y, estimates) {
  n <- length(x)
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Hhat <- empirical_cdf(u, v)
  
  # Compute derivatives at theta=par
  C_vals <- Clayton_Copula(u, v, par)
  dC_vals <- dClayton_Copula_theta(u, v, par)
  d2C_vals <- d2Clayton_Copula_theta(u, v, par)
  
  # One-step formula for WLSE 
  w <- 1 / (Hhat * (1 - Hhat) +  1e-6)   # Variance-stabilizing weights
  B <- 2 * sum( w * (C_vals - Hhat) * dC_vals )
  C_val <- sum( w * (dC_vals^2 + (C_vals - Hhat)*d2C_vals) )
  theta_hat <- par - B / (2*C_val)
  
  return(theta_hat)
}
