#' Fitting SSR Model with MWD Marginals via Clayton Copula
#' @title Estimating parameters of SSR Model with MWD Marginals via Clayton Copula
#' 
#' @description
#' Estimates the parameters of all model parameters such as \eqn{(a_1, b_1,\lambda_1)} for the strength variable \eqn{X},
#' \eqn{(a_2, b_2, \lambda_2)} for the stress variable \eqn{Y}, and \eqn{\theta} for the Clayton copula parameter.
#' 
#' Fits a dependent stress–strength reliability (SSR) model in which both stress
#' and strength marginal distributions follow the Modified Weibull Distribution (MWD),
#' and dependence is modeled using a Clayton copula.
#'
#' The function estimates marginal and copula parameters using several classical
#' methods, including Maximum Likelihood Estimation (MLE), Least Squares Estimation (LSE),
#' Weighted Least Squares Estimation (WLSE), and Maximum Product of Spacings (MPS).
#' 
#' @import stats
#'  
#' @name fit.SSR.ClaytonMWD 
#'  
#' @param data A list with two numeric vectors: \eqn{X} (strength) and \eqn{Y} (stress).
#' @param ACI Logical; if TRUE, \eqn{95\%} asymptotic confidence intervals based on MLEs are computed.
#' @param bootstrap Logical; if TRUE, \eqn{95\%} parametric bootstrap intervals are computed.
#' @param B Number of bootstrap samples.
#' @param seed Integer seed for reproducibility.
#' @param one.step Logical; if TRUE, one-step LSE and WLSE methods are used for estimating \eqn{\theta}.
#' @param alpha Numeric; significance level for confidence intervals (e.g., \eqn{0.05} for \eqn{95\%} CI).
#' @param verbose Logical; if \code{TRUE}, progress and intermediate
#' results from the optimization procedure are printed. Default is \code{FALSE}.
#' 
#' @details
#' Returns point and interval estimates of model parameters using
#' MLE, LSE, WLSE, and MPS methods.
#' Further details can be found in Kızılaslan (2026).
#' 
#' 
#' @return A list containing:
#' \item{all.results}{Point estimates of all model parameters.}
#' \item{theta.Ktau}{Kendall's Tau estimate of \eqn{theta}.}
#' \item{seed}{Used seed in the analysis.}
#' \item{data}{Used data in the analysis.}
#' \item{ACI.parameters}{If ACI=TRUE, it represents the lower and upper bounds of the \eqn{100(1-\alpha)%} 
#' ACI of the parameters with length.}
#' \item{boot.mle}{If bootstrap=TRUE, it represents the lower and upper bounds of the \eqn{100(1-\alpha)%} 
#' parametric bootstrap intervals of 
#' the parameters based on MLE.}
#' \item{boot.lse}{If bootstrap=TRUE, it represents the lower and upper bounds of the \eqn{100(1-\alpha)%} 
#' parametric bootstrap intervals of 
#' the parameters based on LSE.}
#' \item{boot.wlse}{If bootstrap=TRUE, it represents the lower and upper bounds of the \eqn{100(1-\alpha)%} 
#' parametric bootstrap intervals of 
#' the parameters based on WLSE.}
#' \item{boot.mps}{If bootstrap=TRUE, it represents the lower and upper bounds of the \eqn{100(1-\alpha)%} 
#' parametric bootstrap intervals of 
#' the parameters based on MPS.}
#' \item{boot.samples}{A list containing of all the bootstrap samples acrross all the parameters and methods.}
#'
#' @references 
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#'
#' @examples
#' data = list(X = TerkosDam, Y = OmerliDam)
#' fit.SSR = fit.SSR.ClaytonMWD(data,  ACI = TRUE, bootstrap = FALSE, B = 1000,
#'                              seed = 2026, one.step = TRUE, alpha = 0.05)
#' print(fit.SSR)
#' 
#' @export
fit.SSR.ClaytonMWD <- function(data, ACI = FALSE, bootstrap = FALSE, B = NULL, seed = NULL, 
                               one.step = TRUE, alpha = 0.05, verbose = FALSE){

  if(!is.null(seed)) set.seed(seed) else NULL
  n <- length(data$X)
  lower <- c(1e-5, 1e-5, 1e-5) 
  upper <- c(Inf, Inf, Inf)
  init.X <- runif(3) 
  init.Y <- runif(3)  
  
  fit.mleX <- fitMWD(data = data$X, est.method = "mle", opt.method = "L-BFGS-B", starts = init.X, 
                     lower = lower, upper = upper, 
                     hessian = ifelse(ACI==TRUE, TRUE, FALSE) )
  fit.mleY <- fitMWD(data = data$Y, est.method = "mle", opt.method = "L-BFGS-B", starts = init.Y, 
                     lower = lower, upper = upper,
                     hessian = ifelse(ACI==TRUE, TRUE, FALSE) )
  mle.estimates <- as.list(setNames( c(unname(fit.mleX$estimates), unname(fit.mleY$estimates)),
                                     c("a1", "b1", "lambda1", "a2", "b2", "lambda2") ))
  
  fit.lseX <- fitMWD(data = data$X, est.method = "lse", opt.method = "L-BFGS-B", starts = init.X, 
                     lower=lower, upper=upper, hessian = F )
  fit.lseY <- fitMWD(data = data$Y, est.method = "lse", opt.method = "L-BFGS-B", starts = init.Y, 
                     lower=lower, upper=upper, hessian = F )
  lse.estimates <- as.list(setNames( c(unname(fit.lseX$estimates), unname(fit.lseY$estimates)),
                                     c("a1", "b1", "lambda1", "a2", "b2", "lambda2")  ))
  
  fit.wlseX <- fitMWD(data = data$X, est.method = "wlse", opt.method = "L-BFGS-B", starts = init.X, 
                      lower=lower, upper=upper, hessian = F )
  fit.wlseY <- fitMWD(data = data$Y, est.method = "wlse", opt.method = "L-BFGS-B", starts = init.Y, 
                      lower=lower, upper=upper, hessian = F )
  
  wlse.estimates <- as.list(setNames( c(unname(fit.wlseX$estimates), unname(fit.wlseY$estimates)),
                                      c("a1", "b1", "lambda1", "a2", "b2", "lambda2") ))
  
  fit.mpsX <- fitMWD(data = data$X, est.method = "mps", opt.method = "L-BFGS-B", starts = init.X, 
                     lower=lower, upper=upper, hessian = F)
  fit.mpsY <- fitMWD(data = data$Y, est.method = "mps", opt.method = "L-BFGS-B", starts = init.Y, 
                     lower=lower, upper=upper, hessian = F )
  mps.estimates <- as.list(setNames( c(unname(fit.mpsX$estimates), unname(fit.mpsY$estimates)),
                                     c("a1", "b1", "lambda1", "a2", "b2", "lambda2") ))
  
  # Estimates of Clayton dependecy parameter "theta"
  init.theta <- theta_Ktau_estimate(data) # Kendals tau estimate is used as initial 
  fit.theta.mle <- fitClayton(x=data$X, y=data$Y, est.method="mle", opt.method="L-BFGS-B", start=init.theta, 
                              estimates = mle.estimates,
                              lower = 1e-5, upper=Inf, hessian = ifelse(ACI==TRUE, TRUE, FALSE) )

  if(one.step){
    fit.theta.lse <- fit.theta.wlse <- list()
    fit.theta.lse$estimate <- LSE_clayton_onestep(par=init.theta, x=data$X, y=data$Y, estimates = lse.estimates)
    fit.theta.wlse$estimate <- WLSE_clayton_onestep(par=init.theta, x=data$X, y=data$Y, estimates = wlse.estimates)
  }else{
    fit.theta.lse <- fitClayton(x=data$X, y=data$Y, est.method="lse", opt.method="L-BFGS-B", start=init.theta, 
                                estimates = lse.estimates, lower = 1e-5, upper=Inf )
    fit.theta.wlse <- fitClayton(x=data$X, y=data$Y, est.method="wlse", opt.method="L-BFGS-B", start=init.theta, 
                                 estimates = wlse.estimates, lower = 1e-5, upper=Inf )
  }
  fit.theta.mps <- fitClayton(x=data$X, y=data$Y, est.method="mps", opt.method="L-BFGS-B", start=init.theta, estimates = mps.estimates,
                              lower = 1e-5, upper=Inf )
  
  mle.estimates$theta <- fit.theta.mle$estimate
  lse.estimates$theta <- fit.theta.lse$estimate
  wlse.estimates$theta <- fit.theta.wlse$estimate
  mps.estimates$theta <- fit.theta.mps$estimate
  
  Rmle <- Reliability_Clayton_MWD(mle.estimates$a1,mle.estimates$b1, mle.estimates$lambda1, 
                                  mle.estimates$a2, mle.estimates$b2, mle.estimates$lambda2, mle.estimates$theta)$value
  
  Rlse <- Reliability_Clayton_MWD(lse.estimates$a1, lse.estimates$b1, lse.estimates$lambda1, 
                                  lse.estimates$a2, lse.estimates$b2, lse.estimates$lambda2, lse.estimates$theta)$value
  
  Rwlse <- Reliability_Clayton_MWD(wlse.estimates$a1, wlse.estimates$b1, wlse.estimates$lambda1, 
                                   wlse.estimates$a2, wlse.estimates$b2, wlse.estimates$lambda2, wlse.estimates$theta)$value
  
  Rmps <- Reliability_Clayton_MWD(mps.estimates$a1, mps.estimates$b1, mps.estimates$lambda1, 
                                  mps.estimates$a2, mps.estimates$b2, mps.estimates$lambda2, mps.estimates$theta)$value
  
  mle.estimates$R = Rmle
  lse.estimates$R = Rlse
  wlse.estimates$R = Rwlse
  mps.estimates$R = Rmps
  
  # Asymptotic confidence interval based on MLE
  if(ACI){
    parameters.mle <- as.numeric(mle.estimates[1:7])
    hessian.matrix.mle <- matrix( rbind( cbind(fit.mleX$opt.fit$hessian, diag(3)*0 ),
                                         cbind( diag(3)*0, fit.mleY$opt.fit$hessian ) ), nrow = 6, ncol = 6 )
    Sigma_hat <- solve(hessian.matrix.mle)
    se_parameters <- sqrt(diag(Sigma_hat)) # for the marginal parameters
    
    # variance for the copula parameter theta based on Joe (2005)(Journal of Multivariate Analysis 94,401–419).
    M.matrix <- rbind( cbind( hessian.matrix.mle, 0 ), c(rep(0,6), fit.theta.mle$opt.fit$hessian ) )
    
    grad2_l3 <- grad2_nll3(par = parameters.mle, x = data$X, y = data$Y) 
    # 2nd order derivatives of l3(\theta) wrt theta and one of (a1,b1,lambda1,a2,b2,lambda2) at MLEs
    D <- rbind( cbind( hessian.matrix.mle, 0 ), c(grad2_l3)  )
    D_inv <- solve(D)
    V <- D_inv %*% M.matrix %*% t(D_inv)
    se_theta <- sqrt( diag( V ) )[7] 
    #for the copula parameter theta, and others are the same with obtained by marginal Hessian matrices as given in Joe(2005).
    se_parameters <- c(se_parameters, se_theta)
    
    z_alfa2 <- qnorm(1-alpha/2,0,1)
    ACI.parameters <- rbind( lower = pmax(0, parameters.mle - z_alfa2 * se_parameters),
                             upper = parameters.mle + z_alfa2 * se_parameters)

    grad_R <- numDeriv::grad( R_Clayton_MWD, x = parameters.mle )
    var_R <- t(grad_R) %*% V %*% grad_R
    se_R <- sqrt(var_R)
    ACI.R <- c( lower = max(0,Rmle - z_alfa2 * se_R), upper = min(Rmle + z_alfa2 * se_R,1) )
    
    ACI.parameters <- cbind(ACI.parameters, ACI.R)
    ACI.parameters <- rbind(ACI.parameters, ACI.parameters[2,] - ACI.parameters[1,])
    rownames(ACI.parameters)[3] <- "length"
    colnames(ACI.parameters) <- c("a1","b1","lambda1","a2","b2","lambda2","theta","R")
  }
  
  # interval estimates
  if(bootstrap){
    
    # sequential bootstrap for each method
    boot.mle <- parametric_bootstrap(est.method = "mle", opt.method ="L-BFGS-B", boot.estimates = mle.estimates, 
                                     n=n, B = B, seed = seed, one.step = one.step, alpha = alpha)
    boot.mle.ci <- boot.mle$parameters.quantiles

    boot.lse <- parametric_bootstrap(est.method = "lse", opt.method ="L-BFGS-B", boot.estimates = lse.estimates,
                                     n=n, B = B, seed = seed, one.step = one.step, alpha = alpha)
    boot.lse.ci <- boot.lse$parameters.quantiles
    
    boot.wlse <- parametric_bootstrap(est.method = "wlse", opt.method ="L-BFGS-B", boot.estimates = wlse.estimates,
                                      n=n, B = B, seed = seed, one.step = one.step, alpha = alpha)
    boot.wlse.ci <- boot.wlse$parameters.quantiles
    
    boot.mps <- parametric_bootstrap(est.method = "mps", opt.method ="L-BFGS-B", boot.estimates = mps.estimates,
                                     n=n, B = B, seed = seed, one.step = one.step, alpha = alpha)
    boot.mps.ci <- boot.mps$parameters.quantiles

  }
  
  all.results = rbind(mle = as.numeric(mle.estimates),  lse = as.numeric(lse.estimates), 
                      wlse = as.numeric(wlse.estimates), mps = as.numeric(mps.estimates) )
  colnames(all.results) = c("a1", "b1", "lambda1", "a2", "b2", "lambda2", "theta", "R")
  R.results = rbind(RMLE = Rmle, RLSE = Rlse, RWLSE = Rwlse, RMPS = Rmps)

  out <- list(all.results = all.results, theta.Ktau = init.theta, seed = seed, data = data, alpha = alpha )
  class(out) <- "SSRfit"
  
  if(bootstrap){ 
    out$boot.mle = boot.mle.ci
    out$boot.lse = boot.lse.ci
    out$boot.wlse = boot.wlse.ci
    out$boot.mps = boot.mps.ci
    out$boot.samples = list(mle=boot.mle$boot.results, lse=boot.lse$boot.results, wlse=boot.wlse$boot.results, mps=boot.mps$boot.results )
  }
  if(ACI){
    out$ACI.parameters = ACI.parameters
  }
  return(out)
  
}
