#'
#' @title Parametric Bootstrap Confidence Intervals
#' 
#' @description
#' 
#' Computes parametric bootstrap confidence intervals for unknown model parameters
#' and reliability \eqn{R}, based on maximum likelihood estimation (MLE),
#' least squares estimation (LSE), weighted least squares estimation (WLSE),
#' and maximum product of spacing estimation (MPS).
#'
#' @name parametric_bootstrap 
#' 
#' @import stats
#' @importFrom doRNG registerDoRNG
#' 
#' @param est.method Character string specifying the estimation method used.
#'  Options include \code{"MLE"}, \code{"LSE"}, \code{"WLSE"}, and \code{"MPS"}.
#'  
#' @param opt.method Character string specifying the optimization method used in \code{optim}.
#' Common options include \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"},
#' \code{"L-BFGS-B"}, \code{"SANN"}, and \code{"Brent"}.
#'   
#' @param boot.estimates A named list of initial parameter estimates.
#'  The elements \eqn{(a_1, b_1, \lambda_1)} correspond to the strength variable,
#'  \eqn{(a_2, b_2, \lambda_2)} correspond to the stress variable, and
#'  \eqn{\theta} is the Clayton copula dependence parameter.
#'
#' @param n Integer. Sample size.
#'
#' @param B Integer. Number of bootstrap replications.
#'
#' @param seed Integer. Random seed for reproducibility.
#'
#' @param one.step Logical. If \code{TRUE}, one-step LSE and WLSE estimators
#'  are used for \eqn{\theta}.
#'
#' @param alpha Numeric. Significance level for confidence intervals
#'   (e.g., \code{0.05} for a \eqn{95\%} confidence interval).  
#'   
#' @details 
#' This function implements a parametric bootstrap percentile method to construct
#' confidence intervals for unknown parameters and reliability \eqn{R}
#' under different estimation methods (MLE, LSE, WLSE, and MPS).
#' 
#' Further theoretical details are provided in Kizilaslan (2026).
#' 
#' @return A list containing:
#' \item{parameters.quantiles}{A numeric matrix with lower and upper
#'   \eqn{100(1-\alpha)\%} bootstrap percentile confidence limits.}
#' \item{boot.results}{A matrix of bootstrap estimates for all parameters
#'   over \eqn{B} replications.}
#'   
#' 
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @export
parametric_bootstrap <- function(est.method, opt.method, boot.estimates, n, B = 1000, 
                                 seed = NULL, one.step = TRUE, alpha = 0.05){
  
  if(!is.null(seed)) registerDoRNG(seed)
  
  boot.results <- matrix(NA_real_, nrow = B, ncol = length(boot.estimates))
  i <- 1
  
  while(i <= B){
    res <- tryCatch(
      parametric_bootstrap_step(est.method = est.method, opt.method = opt.method, 
                                boot.estimates = boot.estimates, n = n,
                                one.step = one.step),
      error = function(e) NULL
    )
    
    if(!is.null(res) && !any(is.na(res))){
      boot.results[i, ] <- as.numeric(res)
      i <- i + 1  # only move forward if success
    }
    # else retry the same iteration
  }
  
  parameters.quantiles <- apply(boot.results, 2, quantile, c(alpha/2, 1-alpha/2))
  parameters.quantiles <- rbind(parameters.quantiles, parameters.quantiles[2,] - parameters.quantiles[1,] )
  rownames(parameters.quantiles) <- c("lower","upper","length")
  colnames(parameters.quantiles) <- c("a1","b1","lambda1","a2","b2","lambda2","theta","R")
  
  return(list( parameters.quantiles = parameters.quantiles, boot.results = boot.results ) )
}
#'
#' @noRd
parametric_bootstrap_step <-function(est.method, opt.method, boot.estimates, n, seed = NULL, one.step = TRUE){
  
  if (!est.method %in% c("mle", "lse", "wlse", "mps"))
    stop("Estimation method class misspelled. Please check it.")
  
  if (!opt.method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
    stop("Optimization method class misspelled. Please check it.")
  
  if(!is.null(seed)) set.seed(seed) else NULL
  data.bootstrap <- rMweibull_Clayton(n, boot.estimates$a1, boot.estimates$b1, boot.estimates$lambda1, 
                                      boot.estimates$a2, boot.estimates$b2, boot.estimates$lambda2, boot.estimates$theta)
  lower <- c(1e-5, 1e-5, 1e-5);  upper <- c(Inf, Inf, Inf)
  init.X <- runif(1,0.5,1.5)*c(boot.estimates$a1, boot.estimates$b1, boot.estimates$lambda1)
  init.Y <- runif(1,0.5,1.5)*c(boot.estimates$a2, boot.estimates$b2, boot.estimates$lambda2) 

  
  if(est.method == "mle"){
    fit.mleX <- fitMWD(data = data.bootstrap$X, est.method = "mle", opt.method = opt.method, starts = init.X, lower = lower, upper = upper)
    fit.mleY <- fitMWD(data = data.bootstrap$Y, est.method = "mle", opt.method = opt.method, starts = init.Y, lower = lower, upper = upper)
    estimates <- as.list(setNames( c(unname(fit.mleX$estimates), unname(fit.mleY$estimates)),
                                     c("a1", "b1", "lambda1", "a2", "b2", "lambda2") ))
  
    init.theta <- theta_Ktau_estimate(data.bootstrap) # runif(1,0.01,5)
    fit.theta.mle <- fitClayton(x=data.bootstrap$X, y=data.bootstrap$Y, est.method="mle", opt.method=opt.method, start=init.theta, 
                               estimates = estimates, lower = 1e-5, upper=Inf )
    estimates$theta <- fit.theta.mle$estimate
    Rmle <- Reliability_Clayton_MWD(estimates$a1, estimates$b1, estimates$lambda1, 
                                   estimates$a2, estimates$b2, estimates$lambda2, estimates$theta)$value
    estimates$R <- Rmle
  }
  
  if(est.method == "lse"){
    
    fit.lseX <- fitMWD(data = data.bootstrap$X, est.method = "lse", opt.method = opt.method, starts = init.X, lower=lower, upper=upper )
    fit.lseY <- fitMWD(data = data.bootstrap$Y, est.method = "lse", opt.method = opt.method, starts = init.Y, lower=lower, upper=upper )
    estimates <- as.list(setNames( c(unname(fit.lseX$estimates), unname(fit.lseY$estimates)),
                                       c("a1", "b1", "lambda1", "a2", "b2", "lambda2")  ))
    init.theta <- theta_Ktau_estimate(data.bootstrap) 
    if(one.step){
      fit.theta.lse <- list()
      fit.theta.lse$estimate <- LSE_clayton_onestep(par=init.theta, x=data.bootstrap$X, y=data.bootstrap$Y, estimates = estimates)
    }else{
      fit.theta.lse <- fitClayton(x=data.bootstrap$X, y=data.bootstrap$Y, est.method="lse", opt.method=opt.method, start=init.theta, 
                                  estimates = estimates, lower = 1e-5, upper=Inf )
    }
    estimates$theta <- fit.theta.lse$estimate
    Rlse <- Reliability_Clayton_MWD(estimates$a1, estimates$b1, estimates$lambda1, 
                                    estimates$a2, estimates$b2, estimates$lambda2, estimates$theta)$value
    estimates$R <- Rlse
  }
  
  if(est.method == "wlse"){
    
    fit.wlseX <- fitMWD(data = data.bootstrap$X, est.method = "wlse", opt.method = opt.method, starts = init.X, lower=lower, upper=upper )
    fit.wlseY <- fitMWD(data = data.bootstrap$Y, est.method = "wlse", opt.method = opt.method, starts = init.Y, lower=lower, upper=upper )
    estimates <- as.list(setNames( c(unname(fit.wlseX$estimates), unname(fit.wlseY$estimates)),
                                   c("a1", "b1", "lambda1", "a2", "b2", "lambda2")  ))
    init.theta <- theta_Ktau_estimate(data.bootstrap) 
    if(one.step){
      fit.theta.wlse <- list()
      fit.theta.wlse$estimate <- WLSE_clayton_onestep(par=init.theta, x=data.bootstrap$X, y=data.bootstrap$Y, estimates = estimates)
    }else{
      fit.theta.wlse <- fitClayton(x=data.bootstrap$X, y=data.bootstrap$Y, est.method="wlse", opt.method=opt.method, start=init.theta,
                                 estimates = estimates, lower = 1e-5, upper=Inf )
    }
    estimates$theta <- fit.theta.wlse$estimate
    Rwlse <- Reliability_Clayton_MWD(estimates$a1, estimates$b1, estimates$lambda1, 
                                    estimates$a2, estimates$b2, estimates$lambda2, estimates$theta)$value
    estimates$R <- Rwlse
  }
  
  if(est.method == "mps"){
    
    fit.mpsX <- fitMWD(data = data.bootstrap$X, est.method = "mps", opt.method = opt.method, starts = init.X, lower=lower, upper=upper )
    fit.mpsY <- fitMWD(data = data.bootstrap$Y, est.method = "mps", opt.method = opt.method, starts = init.Y, lower=lower, upper=upper )
    estimates <- as.list(setNames( c(unname(fit.mpsX$estimates), unname(fit.mpsY$estimates)),
                                   c("a1", "b1", "lambda1", "a2", "b2", "lambda2")  ))
    init.theta <- theta_Ktau_estimate(data.bootstrap) 
    fit.theta.mps <- fitClayton(x=data.bootstrap$X, y=data.bootstrap$Y, est.method="mps", opt.method=opt.method, start=init.theta,
                                estimates = estimates, lower = 1e-5, upper=Inf )
    estimates$theta <- fit.theta.mps$estimate
    Rmps <- Reliability_Clayton_MWD(estimates$a1, estimates$b1, estimates$lambda1, 
                                     estimates$a2, estimates$b2, estimates$lambda2, estimates$theta)$value
    estimates$R <- Rmps
  }
  
  return(estimates)
}
