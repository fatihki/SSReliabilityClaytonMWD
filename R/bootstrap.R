#'
#' @title Parametric Bootstrap Confidence Intervals
#' @description
#' It gives the parametric bootstrap confidence interval estimation for all the unknown parameters
#' along with the reliability \eqn{R} based on MLE, LSE, WLSE and MPS methods.
#' 
#' @import doRNG stats
#' 
#' @name parametric_bootstrap 
#' 
#' @param est.method Used method for estimating the parameters such as maximum likelihood estimate "MLE",
#' the least square estimation "LSE", the weighted least square estimation "WLSE", and
#' the maximum product of spacing estimates "MPS".
#' @param opt.method The optimization method for \code{optim} function such as "Nelder-Mead", "BFGS", "CG", 
#' "L-BFGS-B", "SANN" and "Brent" to be used for estimating the parameters.
#' @param boot.estimates A list containing the estimates of the unknown parameters.
#' The elements \eqn{(a_1, b_1, \lambda_1)} correspond to the strength variable,
#' \eqn{(a_2, b_2, \lambda_2)} correspond to the stress variable, and \eqn{\theta} corresponds to Clayton copula parameter.
#' @param n sample size
#' @param B Number of bootstrap samples.
#' @param seed Integer seed for reproducibility.
#' @param one.step Logical; if TRUE, one-step LSE and WLSE methods are used for estimating \eqn{\theta}.
#' @param alpha Numeric; significance level for confidence intervals (e.g., \eqn{0.05} for \eqn{95\%} CI).
#' 
#' @details 
#' The parametric bootstrap percentile method is employed to be able to construct bootstrap confidence intervals
#' for unknown parameters and \eqn{R} as well based on MLE, LSE, WLSE and MPSE methods.
#' Further details can be found in Kızılaslan (2026).
#' 
#' @return A list containing:
#' \item{parameters.quantiles}{It gives the \eqn{100(1-\alpha)%} lower and upper bounds.}
#' \item{boot.results}{It represents all the estimates values of the parameters through B bootstrap samples.}
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
    if ( is.null(fit.mleX)) {
    next
      }
  
   fit.mleY <- fitMWD(data = data.bootstrap$Y, est.method = "mle", opt.method = opt.method, starts = init.Y, lower = lower, upper = upper)
   if ( is.null(fit.mleY)) {
    next
   }
   estimates <- as.list(setNames( c(unname(fit.mleX$estimates), unname(fit.mleY$estimates)),
                                     c("a1", "b1", "lambda1", "a2", "b2", "lambda2") ))
  
   init.theta <- theta_Ktau_estimate(data.bootstrap) # runif(1,0.01,5)
   fit.theta.mle <- fitClayton(x=data.bootstrap$X, y=data.bootstrap$Y, est.method="mle", opt.method=opt.method, start=init.theta, 
                               estimates = estimates, lower = 1e-5, upper=Inf )
   if ( is.null(fit.theta.mle)) {
     next
   }
   estimates$theta <- fit.theta.mle$estimate
   Rmle <- Reliability_Clayton_MWD(estimates$a1, estimates$b1, estimates$lambda1, 
                                   estimates$a2, estimates$b2, estimates$lambda2, estimates$theta)$value
   estimates$R = Rmle
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
    estimates$R = Rlse
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
    estimates$R = Rwlse
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
    estimates$R = Rmps
  }
  
  return(estimates)
}
