#'  
#' @title Estimating parameter of the Clayton copula
#' 
#' @description
#' Estimates the dependence parameter \eqn{\theta} of the Clayton copula.
#' 
#' @name fitClayton
#'
#' @param x Vector of observations for the strength variable \eqn{X}.
#' @param y Vector of observations for the stress variable \eqn{Y}. 
#' @param est.method Used method for estimating the parameters such as maximum likelihood estimate "MLE",
#' the least square estimation "LSE", the weighted least square estimation "WLSE", and
#' the maximum product of spacing estimates "MPS".
#' @param opt.method The optimization method for \code{optim} function such as "Nelder-Mead", "BFGS", "CG", 
#' "L-BFGS-B", "SANN" and "Brent" to be used for estimating the parameters.
#' @param start Initial value for \eqn{\theta} parameter to be optimized over.
#' @param estimates A list containing estimates of all model paramters \eqn{a_1,b_1,\lambda_1, a_2, b_2, \lambda_2}.
#' @param lower Numeric vector specifying the lower bounds of the parameters
#' for bounded optimization methods (e.g., \code{"L-BFGS-B"}). Ignored if \code{NULL}.
#' @param upper Numeric vector specifying the upper bounds of the parameters
#' for bounded optimization methods. Ignored if \code{NULL}.
#' @param verbose Logical; if \code{TRUE}, progress and intermediate
#' results from the optimization procedure are printed. Default is \code{FALSE}.
#' @param ... Further arguments to be passed to \code{optim} function such as lower and upper limits, hessian, etc.
#' 
#' @details 
#' The joint distribution function of the two-dimensional Clayton copula is
#' \deqn{
#' C(u,v;\theta) = \left(u^{-\theta} + v^{-\theta} - 1\right)^{-1/\theta},
#' }
#' where the dependence parameter satisfies \eqn{\theta > 0}.
#' It is estimated using several classical methods:
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
#' Further details can be found in Kızılaslan (2026).
#' 
#' @return A list containing:
#' \item{estimate}{Estimate of the Clayton copula parameter, \eqn{\theta}.}
#' \item{opt.fit}{Full optimization result.}
#' 
#' @references
#' Kızılaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' arXiv preprint. Available at
#' \href{https://arxiv.org/abs/????}{https://arxiv.org/abs/????}.
#' 
#' @export
#' 
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
      message("Optimization failed")
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
        message("Optimization failed")
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
        message("Optimization failed")
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
        message("Optimization failed")
        NULL
      }
      )

  }
  
  if (is.null(out)) {
    message("Optimization failed — exiting this run.")
    return(NULL)  # or stop() if you want to terminate entirely
  }
  
  par.est <- out$par
  n <- length(data)
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
#' since the original obj function is maximized, it is multiplying with -1 for minimizing
#' @keywords internal
mps_clayton <- function(par, x, y, estimates) {
  if(par < 0 ) return(-Inf)  
  n <- length(x) 
  u <- pMweibull(x, estimates$a1, estimates$b1, estimates$lambda1)
  v <- pMweibull(y, estimates$a2, estimates$b2, estimates$lambda2)
  Ccopula <- Clayton_Copula(u, v, par)
  Ccopula <- sort(Ccopula)
  D3 <- diff( c(0, Ccopula, 1) )

  return( -sum(log(D3))/(n+1) )
}
#' Kendall's Tau Estimate of Clayton copula parameter theta
#' @name theta_Ktau_estimate
#' 
#' @title Kendall's Tau Estimate of Clayton copula parameter
#' 
#' @param data A list with two numeric vectors: \eqn{X} (strength) and \eqn{Y} (stress).
#'
#' @return Moment estimate of \eqn{\theta} parameter based on the Kendall's \eqn{\tau}.
#' 
#' @examples
#' set.seed(123)
#' n <- 50
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 5 # 1, 2, 3, 4, 5
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
#'
# -------------------------------
# One-step LSE estimate
# -------------------------------
#' @name LSE_clayton_onestep
#' 
#' @title One-step LSE estimate of the Clayton copula parameter
#' 
#' @param par Numeric value of the dependence parameter \eqn{\theta} of the Clayton copula.
#' @param x Vector of observations for the strength variable \eqn{X}.
#' @param y Vector of observations for the stress variable \eqn{Y}.
#' @param estimates A list containing estimates of the model parameters 
#' \eqn{a_1, b_1, \lambda_1, a_2, b_2, \lambda_2}.
#' 
#' @details
#' Further details are provided in Kızılaslan (2026).
#' 
#' @return The one-step LSE of \eqn{\theta} parameter.
#' 
#' @references
#' Kızılaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/????}{arXiv:????}
#' 
#' @export
LSE_clayton_onestep <- function( par, x, y, estimates) {
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
#' @name WLSE_clayton_onestep
#' 
#' @title One-step WLSE estimate of the Clayton copula parameter
#' 
#' @param par Numeric value of the dependence parameter \eqn{\theta} of the Clayton copula.
#' @param x Vector of observations for the strength variable \eqn{X}.
#' @param y Vector of observations for the stress variable \eqn{Y}.
#' @param estimates A list containing estimates of the model parameters 
#' \eqn{a_1, b_1, \lambda_1, a_2, b_2, \lambda_2}.
#' 
#' @details
#' Further details are provided in Kızılaslan (2026).
#' 
#' @return The one-step WLSE of \eqn{\theta} parameter.
#' 
#' @references
#' Kızılaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/????}{arXiv:????}
#' 
#' @export
#' 
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
